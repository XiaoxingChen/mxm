#if !defined(__SPATIAL_BVH_H__)
#define __SPATIAL_BVH_H__

#include "spatial_aabb.h"
#include "geometry_primitive.h"
#include "spatial_tree_base.h"
#include <stack>
#include <memory>
#include <map>


namespace mxm
{

namespace bvh
{

enum HitType
{
    eClosestHit = 0,
    eMultiHit = 1,
    eAnyHit = 2,
    eHitTypeNum
};

template<typename DType=float>
struct Node
{
    Node<DType>(size_t dim): aabb(dim), is_leaf(false){}
    AABB<DType> aabb;
    std::vector<size_t> children_index_buffer;
    std::vector<size_t> primitive_index_buffer;
    bool is_leaf;
    size_t prim_num;
};

template<typename DType=float>
class BaseTree
{
public:
    virtual size_t dim() const = 0;
    virtual Matrix<DType> primitive(size_t idx) const = 0;
    virtual size_t primitiveSize() const = 0;
    virtual void build(size_t primitive_per_leaf=4, bool verbose=true);
    const std::vector<bvh::Node<DType>>& nodeBuffer() const { return node_buffer_;}

protected:
    std::vector<bvh::Node<DType>> node_buffer_;
};

class PrimitiveMeshTree: public BaseTree<>
{
public:
    struct HitRecord
    {
        Ray<> ray;
        FloatType t;
        size_t prim_idx;
        Vector<FloatType> coeff;
    };

    // virtual void build(size_t primitive_per_leaf=4, bool verbose=true) override;
    virtual size_t dim() const override { return vertex_buffer_->shape(0); }

    PrimitiveMeshTree(
        std::shared_ptr<Mat>& vertex_buffer,
        std::shared_ptr<Matrix<size_t>>& vertex_index_buffer)
        :vertex_buffer_(vertex_buffer), vertex_index_buffer_(vertex_index_buffer)
        {}

    PrimitiveMeshTree()
        :vertex_buffer_(nullptr), vertex_index_buffer_(nullptr)
        {}

    virtual Mat primitive(size_t idx) const override { return getPrimitive(*vertex_buffer_, *vertex_index_buffer_, idx); }

    size_t multiHit(const Ray<>& ray) const;
    std::vector<HitRecord> hit(const Ray<>& ray, HitType hit_type) const;

    // const AABB& aabb() const { return node_buffer_.at(0).aabb; }
    const Mat & vertexBuffer() const { return *vertex_buffer_; }
    virtual size_t primitiveSize() const override { return vertex_index_buffer_->shape(1); }
    // const std::vector<bvh::Node>& nodeBuffer() const { return node_buffer_;}

private:

    std::shared_ptr<Mat> vertex_buffer_;
    std::shared_ptr<Matrix<size_t>> vertex_index_buffer_;
    // std::vector<bvh::Node> node_buffer_;

};

size_t rayCast(
    const Ray<>& ray,
    const PrimitiveMeshTree& tree,
    HitType hit_type,
    Matrix<float>* p_coeff=nullptr,
    Vector<size_t>* p_prim_idx=nullptr,
    Vector<float>* p_hit_t=nullptr);

class PointCloudTree: public BaseTree<>
{
public:
    PointCloudTree(
        std::shared_ptr<Mat>& point_buffer)
        :point_buffer_(point_buffer) {}
    virtual size_t dim() const override {return point_buffer_->shape(0);}
    size_t pointNum() const {return point_buffer_->shape(1);}

    // void build(size_t point_per_leaf=4, bool verbose=true);
    // const std::vector<bvh::Node>& nodeBuffer() const { return node_buffer_; }
    virtual Mat primitive(size_t idx) const override { return (*point_buffer_)(Col(idx)); }
    virtual size_t primitiveSize() const override { return point_buffer_->shape(1); }
    std::multimap<FloatType, size_t> radiusSearch(const Vec& pt, FloatType radius) const;
    std::multimap<FloatType, size_t> nearestNeighborSearch(const Vec& pt, size_t k) const;

private:
    std::shared_ptr<Mat> point_buffer_;
    // std::vector<bvh::Node> node_buffer_;
};

size_t treeNodeRequirement(size_t leaf_num, size_t child_num);
bool verifyTree(const std::vector<bvh::Node<>>& node_buffer, size_t node_idx);

template<typename DType>
void BaseTree<DType>::build(size_t primitive_per_leaf, bool verbose)
{
    std::vector<size_t> primitive_index_buffer(primitiveSize());
    for(size_t i = 0; i < primitive_index_buffer.size(); i++) primitive_index_buffer.at(i) = i;

    {
        size_t require = treeNodeRequirement(primitive_index_buffer.size() / primitive_per_leaf + 1, 2);
        node_buffer_.clear();
        node_buffer_.reserve(require);
    }

    std::stack<RangeNode> stk;

    stk.push({node_buffer_.size(), {0, primitive_index_buffer.size()}});
    node_buffer_.push_back(Node(dim()));

    auto t_start = std::chrono::system_clock::now();
    while (!stk.empty())
    {
        auto target = stk.top();
        stk.pop();

        auto& target_node = node_buffer_.at(target.node_idx);
        target_node.prim_num = target.range[1] - target.range[0];

        //
        // eliminate leaf node
        if(target.range[1] - target.range[0] <= primitive_per_leaf)
        {
            for(size_t sorted_idx = target.range[0]; sorted_idx < target.range[1]; sorted_idx++)
            {
                auto prim_idx = primitive_index_buffer.at(sorted_idx);
                target_node.primitive_index_buffer.push_back(prim_idx);
                target_node.aabb.extend( primitive(prim_idx) );
            }

            target_node.is_leaf = true;
            continue;
        }

        //
        // deal with internal node
        //

        //
        // calculate AABB
        for(size_t i = target.range[0]; i < target.range[1]; i++)
        {
            const size_t& prim_idx = primitive_index_buffer.at(i);
            target_node.aabb.extend(primitive(prim_idx));
        }

        //
        // find the longest axis
        size_t target_axis = argMax(target_node.aabb.max() - target_node.aabb.min())[0];
        if(verbose)
            std::cout << "range: [" << target.range[0] << "," << target.range[1] << "), axis: " << target_axis << std::endl;

        //
        // sort along the target axis
        size_t mid = (target.range[1] + target.range[0]) / 2;

        std::nth_element(
            primitive_index_buffer.begin() + target.range[0],
            primitive_index_buffer.begin() + mid,
            primitive_index_buffer.begin() + target.range[1],
            [&](const size_t& prim_idx1, const size_t& prim_idx2)
            {
                Mat prim1 = primitive(prim_idx1);
                Mat prim2 = primitive(prim_idx2);

                return mxm::sum(prim1(Row(target_axis))) < mxm::sum(prim2(Row(target_axis)));
            }
        );

        // create children

        stk.push({node_buffer_.size(), {target.range[0], mid}});
        target_node.children_index_buffer.push_back(node_buffer_.size());
        node_buffer_.push_back(Node(dim()));

        stk.push({node_buffer_.size(), {mid, target.range[1]}});
        target_node.children_index_buffer.push_back(node_buffer_.size());
        node_buffer_.push_back(Node(dim()));

    }
    auto t_end = std::chrono::system_clock::now();
    std::cout << "build tree t_cost(s): " << std::chrono::duration<double>(t_end - t_start).count() << std::endl;
    if(!verifyTree(node_buffer_, 0))
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}


} // namespace bvh


} // namespace mxm

#ifdef MXM_HEADER_ONLY
#include "spatial_bvh_inl.h"
#endif



#endif // __SPATIAL_BVH_H__
