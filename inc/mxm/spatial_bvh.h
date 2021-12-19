#if !defined(__SPATIAL_BVH_H__)
#define __SPATIAL_BVH_H__

#include "spatial_aabb.h"
#include "geometry_primitive.h"
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

struct Node
{
    Node(size_t dim): aabb(dim), is_leaf(false){}
    AABB aabb;
    std::vector<size_t> children_index_buffer;
    std::vector<size_t> primitive_index_buffer;
    bool is_leaf;
    size_t prim_num;
};

class BaseTree
{
public:
    virtual size_t dim() const = 0;
    virtual Mat primitive(size_t idx) const = 0;
    virtual size_t primitiveSize() const = 0;
    virtual void build(size_t primitive_per_leaf=4, bool verbose=true);
    const std::vector<bvh::Node>& nodeBuffer() const { return node_buffer_;}

protected:
    std::vector<bvh::Node> node_buffer_;
};

class PrimitiveMeshTree: public BaseTree
{
public:
    struct HitRecord
    {
        Ray ray;
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

    virtual Mat primitive(size_t idx) const override { return getPrimitive(*vertex_buffer_, *vertex_index_buffer_, idx); }

    size_t multiHit(const Ray& ray) const;
    std::vector<HitRecord> hit(const Ray& ray, HitType hit_type) const;

    // const AABB& aabb() const { return node_buffer_.at(0).aabb; }
    const Mat & vertexBuffer() const { return *vertex_buffer_; }
    virtual size_t primitiveSize() const override { return vertex_index_buffer_->shape(1); }
    // const std::vector<bvh::Node>& nodeBuffer() const { return node_buffer_;}

private:

    std::shared_ptr<Mat> vertex_buffer_;
    std::shared_ptr<Matrix<size_t>> vertex_index_buffer_;
    // std::vector<bvh::Node> node_buffer_;

};

class PointCloudTree: public BaseTree
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


} // namespace bvh


} // namespace mxm

#ifdef MXM_HEADER_ONLY
#include "spatial_bvh_inl.h"
#endif



#endif // __SPATIAL_BVH_H__
