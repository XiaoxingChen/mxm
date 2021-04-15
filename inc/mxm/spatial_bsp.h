#if !defined(__SPATIAL_BSP_H__)
#define __SPATIAL_BSP_H__

#include "linalg.h"
#include <memory>
#include <stack>

#include "spatial_tree_base.h"

namespace mxm
{

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

private:
    std::shared_ptr<Mat> point_buffer_;
    std::vector<bsp::Node> node_buffer_;
};

inline void divide(
    std::vector<size_t>::iterator begin,
    std::vector<size_t>::iterator end,
    const Mat& points,
    size_t axis,
    FloatType compare_val)
{
    auto i = begin;
    auto j = std::prev(end, 1);
    while(i < j)
    {
        while(points(axis, *i) <= compare_val) i++;
        while(points(axis, *j) > compare_val) j++;
        if (i >= j) break;
        auto temp = *i;
        *i = *j;
        *j = temp;
    }
}
// todo test divide

template<typename DType>
Vector<DType> binaryToVector(size_t dim, uint32_t bin)
{
    Vector<DType> ret(Vector<DType>::zeros(dim));
    size_t axis = 0;
    while(bin > 0)
    {
        if((bin & 1u) > 0) ret(axis) = DType(1);
        bin >>= 1;
    }
    return ret;
}

inline void PointCloudTree::build(size_t point_per_leaf, bool verbose)
{
    size_t total_node_num = (pointNum() / point_per_leaf + 1) * 2;
    node_buffer_.reserve(total_node_num);
    std::vector<size_t> index_buffer;

    std::stack<RangeNode> stk;

    stk.push({node_buffer_.size(), {0, index_buffer.size()}});

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
                target_node.point_index_buffer.push_back(idx);
                target_node.is_leaf = true;
            }

        }else // internal node
        {
            std::vector<size_t> mid_list(dim() + 1);
            // todo: assign value to mid_list;

            for(size_t i = 0; i < (1ul << dim()); i++)
            {
                stk.push({node_buffer_.size(), {mid_list.at(i), mid_list.at(i+1)}});
                target_node.children_index_buffer.push_back(node_buffer_.size());
                node_buffer_.push_back(Node(dim()));
                node_buffer_.back().min = target_node.min + binaryToVector<FloatType>(dim(), i) * target_node.width * 0.5;
                node_buffer_.back().width = target_node.width * 0.5;
            }

        }
        // for(size_t axis = 0; axis < dim(); axis++)
        // {
        //     divide(
        //         index_buffer.begin() + target.range[0],
        //         index_buffer.begin() + target.range[1],
        //         point_buffer_, axis,

        // }
    }


}

} // namespace bsp


} // namespace mxm



#endif // __SPATIAL_BSP_H__
