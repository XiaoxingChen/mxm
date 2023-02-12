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
    size_t height = static_cast<size_t>(1. + std::log(leaf_num) / std::log(child_num)) + 1;
    if(child_num != 2)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    return (1 << height);
}

MXM_INLINE bool verifyTree(const std::vector<bvh::Node<>>& node_buffer, size_t node_idx)
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

MXM_INLINE size_t PrimitiveMeshTree::multiHit(const Ray<>& ray) const
{
    auto records = hit(ray, eMultiHit);
    return records.size();
}

MXM_INLINE size_t rayCast(
    const Ray<>& ray_in,
    const PrimitiveMeshTree& tree,
    HitType hit_type,
    Matrix<float>* p_coeff,
    Vector<size_t>* p_prim_idx,
    Vector<float>* p_hit_t)
{
    Ray<> ray(ray_in);
    std::stack<size_t> node_idx_stk;
    std::vector<float> data_hit_t;
    // std::vector<float> data_vertex_coeff;
    // Matrix<float> bary_coords({ray.origin().size(), 0});
    Matrix<float> bary_coords({ray.origin().size(), 0});
    std::vector<size_t> data_prim_idx;

    if(!tree.nodeBuffer().at(0).aabb.hit(ray))
        return 0;
    node_idx_stk.push(0);

    while(! node_idx_stk.empty())
    {
        const Node<>& target_node = tree.nodeBuffer().at(node_idx_stk.top());
        node_idx_stk.pop();

        //
        // deal with leaf node
        if(target_node.is_leaf)
        {
            for(const auto& prim_idx: target_node.primitive_index_buffer)
            {
                const auto [hit_t, bary_coord] = rayPrimitiveIntersection(ray.origin(), ray.direction(), tree.primitive(prim_idx));
                if(min(bary_coord) < -eps<float>() || !ray.valid(hit_t))
                    continue;

                if(hit_type == eAnyHit
                || hit_type == eMultiHit
                || (hit_type == eClosestHit && data_hit_t.empty()))
                {
                    data_hit_t.push_back(hit_t);
                    data_prim_idx.push_back(prim_idx);
                    bary_coords = hstack(bary_coords, bary_coord);
                }else if (hit_type == eClosestHit && !data_hit_t.empty() && hit_t < data_hit_t.at(0))
                {
                    data_hit_t.at(0) = hit_t;
                    data_prim_idx.at(0) = prim_idx;
                    bary_coords = hstack(bary_coords, bary_coord);
                }else{ assert(false); }


                if(hit_type == eAnyHit) break;
                if(hit_type == eClosestHit)
                {
                    ray.tMax() = hit_t;
                }

            }
            continue;
        }

        //
        // deal with internal node
        for(const auto& child_idx: target_node.children_index_buffer)
        {
            if(tree.nodeBuffer().at(child_idx).aabb.hit(ray))
                node_idx_stk.push(child_idx);
        }

    }
    size_t hit_cnt = data_hit_t.size();

    if(p_coeff) (*p_coeff) = bary_coords;
    if(p_prim_idx) (*p_prim_idx) = Vector<size_t>(std::move(data_prim_idx));
    if(p_hit_t) (*p_hit_t) = Vector<float>(std::move(data_hit_t));

    return hit_cnt;
}

MXM_INLINE std::vector<PrimitiveMeshTree::HitRecord> PrimitiveMeshTree::hit(const Ray<>& ray_in, HitType hit_type) const
{
    std::vector<HitRecord> records;
    Matrix<float> coeff;
    Vector<size_t> prim_idx;
    Vector<float> hit_t;
    rayCast(ray_in, *this, hit_type, &coeff, &prim_idx, &hit_t);
    for(size_t i = 0; i < prim_idx.size(); i++)
    {
        HitRecord record;
        record.ray = ray_in;
        record.prim_idx = prim_idx(i);
        record.coeff = coeff(Col(i));
        record.t = hit_t(i);
        records.push_back(record);
    }

    return records;
}

MXM_INLINE std::multimap<FloatType, size_t> PointCloudTree::radiusSearch(const Vec& pt, FloatType radius) const
{
    std::stack<Node<>> stk;
    stk.push(nodeBuffer().at(0));

    std::multimap<FloatType, size_t> ret;

    while(!stk.empty())
    {

        auto target_node = stk.top();
        stk.pop();

        // leaf node
        if(target_node.is_leaf)
        {
            for(auto idx : target_node.primitive_index_buffer)
            {
                auto dist = ((*point_buffer_)(Col(idx)) - pt).norm();
                if(dist < radius)
                    ret.insert({dist, idx});
            }
        }else // internal node
        {
            for(auto idx: target_node.children_index_buffer)
            {
                if(distance(nodeBuffer().at(idx).aabb, pt)(0,0) < radius)
                    stk.push(nodeBuffer().at(idx));
            }
        }
    }
    return ret;
}

MXM_INLINE std::multimap<FloatType, size_t> PointCloudTree::nearestNeighborSearch(const Vec& pt, size_t k) const
{
    std::multimap<FloatType, size_t> ret;

    std::stack<size_t> stk;
    stk.push(0);
    while(!stk.empty())
    {
        auto target = nodeBuffer().at(stk.top());
        stk.pop();

        auto curr_max_dist = ret.empty() ? INFINITY : std::prev(ret.end())->first;

        if(target.is_leaf)
        {
            for(auto idx: target.primitive_index_buffer)
            {
                auto dist = ((*point_buffer_)(Col(idx)) - pt).norm();
                if(dist < curr_max_dist)
                    ret.insert({dist, idx});
            }
        }else
        {
            for(auto idx: target.children_index_buffer)
            {
                auto dist = distance(nodeBuffer().at(idx).aabb, pt)(0,0);
                if(dist < curr_max_dist) stk.push(idx);
            }
        }


    while(ret.size() > k)
    {
        ret.erase(std::prev(ret.end()));
    }

    }
    return ret;
}

} // namespace bvh

} // namespace mxm




#endif // __SPATIAL_BVH_INL_H__
