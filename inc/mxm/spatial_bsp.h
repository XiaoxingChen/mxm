#if !defined(__SPATIAL_BSP_H__)
#define __SPATIAL_BSP_H__

#include "linalg.h"
#include "spatial_aabb.h"
#include <memory>
#include <stack>
#include <set>

#include "spatial_tree_base.h"

namespace mxm
{

template<typename DType>
Vector<DType> binaryToVector(size_t dim, uint32_t bin)
{
    Vector<DType> ret(Vector<DType>::zeros(dim));
    size_t axis = 0;
    while(axis < dim)
    {
        if(((bin >> axis) & 1u) > 0) ret(axis) = DType(1);
        axis++;
    }
    return ret;
}

// Generalization of Binary Space Partition:
// Binary tree, quadtree, octree.
namespace bsp
{

struct Node
{
    Node(size_t dim): min(dim), is_leaf(false){}
    Vec min;
    FloatType width;

    std::vector<size_t> children_index_buffer;
    std::vector<size_t> point_index_buffer;
    bool is_leaf;
};

class PointCloudTree
{
public:
    PointCloudTree(
        std::shared_ptr<Mat>& point_buffer)
        :point_buffer_(point_buffer) {}
    size_t dim() const {return point_buffer_->shape(0);}
    size_t pointNum() const {return point_buffer_->shape(1);}

    void build(size_t point_per_leaf=4, bool verbose=true);
    const std::vector<bsp::Node>& nodeBuffer() const { return node_buffer_; }

private:
    std::shared_ptr<Mat> point_buffer_;
    std::vector<bsp::Node> node_buffer_;
};

inline size_t axisPartition(
    std::vector<size_t>& point_index_buffer,
    std::array<size_t, 2> range,
    const Mat& points,
    size_t axis,
    const Vec& thresholds)
{
    auto i = range[0];
    auto j = range[1] - 1;
    while(i < j)
    {
        while(points(axis, point_index_buffer.at(i)) <= thresholds(axis)) i++;
        while(points(axis, point_index_buffer.at(j)) > thresholds(axis)) j--;
        if (i >= j) break;
        auto temp = point_index_buffer.at(i);
        point_index_buffer.at(i) = point_index_buffer.at(j);
        point_index_buffer.at(j) = temp;
    }
    return i;
}


inline std::vector<size_t> fullAxesPartition(
    std::vector<size_t>& point_index_buffer,
    std::array<size_t, 2> range,
    const Mat& points,
    const Vec& thresholds)
{
    size_t n_axes = points.shape(0);
    std::vector<size_t> partitioners{range[0], range[1]};

    struct Context
    {
        std::array<size_t, 2> range;
        size_t axis;
    };
    std::stack<Context> stk;
    stk.push(Context{{range[0], range[1]}, 0});

    // pre-order traversal
    while(!stk.empty())
    {
        auto context = stk.top();
        stk.pop();

        auto partiton_idx = axisPartition(point_index_buffer, context.range, points, context.axis, thresholds);

        partitioners.push_back(partiton_idx);
        if(context.axis + 1 < n_axes)
        {
            stk.push(Context{{partiton_idx, context.range[1]}, context.axis + 1});
            stk.push(Context{{context.range[0], partiton_idx}, context.axis + 1});
        }
    }

    std::sort(partitioners.begin(), partitioners.end());
    return partitioners;
}


inline void PointCloudTree::build(size_t point_per_leaf, bool verbose)
{
    size_t total_node_num = (pointNum() / point_per_leaf + 1) * 2;
    node_buffer_.reserve(total_node_num);
    std::vector<size_t> index_buffer;
    for(size_t i = 0; i < point_buffer_->shape(1); i++) index_buffer.push_back(i);

    std::stack<RangeNode> stk;

    stk.push({node_buffer_.size(), {0, index_buffer.size()}});
    node_buffer_.push_back(Node(dim()));

    {
        AABB aabb(dim());
        aabb.extend(*point_buffer_);
        node_buffer_.back().min = aabb.min();
        node_buffer_.back().width = mxm::min(aabb.max() - aabb.min());
    }


    while(!stk.empty())
    {
        auto target = stk.top();
        stk.pop();

        auto& target_node = node_buffer_.at(target.node_idx);

        // leaf node
        if(target.range[1] - target.range[0] <= point_per_leaf)
        {
            for(size_t idx = target.range[0]; idx < target.range[1]; idx++)
            {
                target_node.point_index_buffer.push_back( index_buffer.at(idx) );
                target_node.is_leaf = true;
            }

        }else // internal node
        {

            Vec thresholds = target_node.min + target_node.width * 0.5;
            auto mid_list = fullAxesPartition(index_buffer, target.range, *point_buffer_, thresholds);

            for(size_t i = 0; i < (1ul << dim()); i++)
            {
                stk.push({node_buffer_.size(), {mid_list.at(i), mid_list.at(i+1)}});
                target_node.children_index_buffer.push_back(node_buffer_.size());
                node_buffer_.push_back(Node(dim()));
                node_buffer_.back().min = target_node.min + binaryToVector<FloatType>(dim(), i) * target_node.width * 0.5;
                node_buffer_.back().width = target_node.width * 0.5;
            }

        }
    }


}

} // namespace bsp


} // namespace mxm



#endif // __SPATIAL_BSP_H__
