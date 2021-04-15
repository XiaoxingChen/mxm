#if !defined(__SPATIAL_BVH_INL_H__)
#define __SPATIAL_BVH_INL_H__

#if defined(MXM_COMPILED_LIB)
#include "spatial_bvh.h"
#endif // MXM_COMPILED_LIB

#include <array>
#include <cmath>
#include <chrono>

#include "spatial_tree_base.h"

namespace mxm
{
namespace bvh
{


// struct RangeNode
// {
//     size_t node_idx;
//     std::array<size_t, 2> range;
// };

MXM_INLINE size_t treeNodeRequirement(size_t leaf_num, size_t child_num=2)
{
    size_t height = static_cast<size_t>(1. + log(leaf_num) / log(child_num)) + 1;
    if(child_num != 2)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    return (1 << height);
}

MXM_INLINE bool verifyTree(const std::vector<bvh::Node>& node_buffer, size_t node_idx)
{
    if(node_buffer.at(node_idx).is_leaf)
        return true;
    for(auto & child_idx: node_buffer.at(node_idx).children_index_buffer)
    {
        if(!node_buffer.at(child_idx).aabb.in(node_buffer.at(node_idx).aabb))
            return false;
        if(!verifyTree(node_buffer, child_idx))
            return false;
    }
    return true;
}

MXM_INLINE void PrimitiveMeshTree::build(size_t primitive_per_leaf, bool verbose)
{
    std::vector<size_t> primitive_index_buffer(vertex_index_buffer_->shape(1));
    for(size_t i = 0; i < primitive_index_buffer.size(); i++) primitive_index_buffer.at(i) = i;

    {
        size_t require = treeNodeRequirement(primitive_index_buffer.size() / primitive_per_leaf + 1);
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
        size_t target_axis = argMax(target_node.aabb.max() - target_node.aabb.min());
        if(verbose)
            std::cout << "range: [" << target.range[0] << "," << target.range[1] << "), axis: " << target_axis << std::endl;

        //
        // sort along the target axis
        size_t mid = (target.range[1] + target.range[0]) / 2;

        #if 0
        std::sort(
            primitive_index_buffer.begin() + target.range[0],
            primitive_index_buffer.begin() + target.range[1],
            [&](const size_t& prim_idx1, const size_t& prim_idx2)
            {
                Mat prim1 = primitive(prim_idx1);
                Mat prim2 = primitive(prim_idx2);

                return prim1(Row(target_axis)).asVector().sum() < prim2(Row(target_axis)).asVector().sum();
            });

        #else
        std::nth_element(
            primitive_index_buffer.begin() + target.range[0],
            primitive_index_buffer.begin() + mid,
            primitive_index_buffer.begin() + target.range[1],
            [&](const size_t& prim_idx1, const size_t& prim_idx2)
            {
                Mat prim1 = primitive(prim_idx1);
                Mat prim2 = primitive(prim_idx2);

                return prim1(Row(target_axis)).asVector().sum() < prim2(Row(target_axis)).asVector().sum();
            }
        );
        #endif

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

MXM_INLINE size_t PrimitiveMeshTree::multiHit(const Ray& ray) const
{
    auto records = hit(ray, eMultiHit);
    return records.size();
}
MXM_INLINE std::vector<PrimitiveMeshTree::HitRecord> PrimitiveMeshTree::hit(const Ray& ray_in, HitType hit_type) const
{
    std::vector<HitRecord> records;
    Ray ray(ray_in);
    std::stack<size_t> node_idx_stk;

    if(!node_buffer_.at(0).aabb.hit(ray))
        return records;
    node_idx_stk.push(0);

    while(! node_idx_stk.empty())
    {
        const Node& target_node = node_buffer_.at(node_idx_stk.top());
        node_idx_stk.pop();

        //
        // deal with leaf node
        if(target_node.is_leaf)
        {
            for(const auto& prim_idx: target_node.primitive_index_buffer)
            {
                auto result = intersectEquation( primitive(prim_idx), ray);
                auto hit_t = result(0);
                if(!validIntersect (result) || !ray.valid(hit_t))
                    continue;

                HitRecord record;
                record.t = hit_t;
                record.prim_idx = prim_idx;
                record.ray = ray_in;


                if(hit_type == eAnyHit)
                {
                    records.push_back(record);
                    return records;
                }
                if(hit_type == eClosestHit)
                {
                    if(records.empty()) records.push_back(record);
                    if(records.front().t > record.t) records.front() = record;
                    ray.tMax() = hit_t;
                }else if(hit_type == eMultiHit)
                {
                    records.push_back(record);
                }

            }
            continue;
        }

        //
        // deal with internal node
        for(const auto& child_idx: target_node.children_index_buffer)
        {
            if(node_buffer_.at(child_idx).aabb.hit(ray))
                node_idx_stk.push(child_idx);
        }

    }

    return records;
}

} // namespace bvh

} // namespace mxm




#endif // __SPATIAL_BVH_INL_H__
