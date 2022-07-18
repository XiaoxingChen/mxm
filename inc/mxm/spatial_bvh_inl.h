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

MXM_INLINE void BaseTree::build(size_t primitive_per_leaf, bool verbose)
{
    std::vector<size_t> primitive_index_buffer(primitiveSize());
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
    std::vector<float> data_vertex_coeff;
    std::vector<size_t> data_prim_idx;

    if(!tree.nodeBuffer().at(0).aabb.hit(ray))
        return 0;
    node_idx_stk.push(0);

    while(! node_idx_stk.empty())
    {
        const Node& target_node = tree.nodeBuffer().at(node_idx_stk.top());
        node_idx_stk.pop();

        //
        // deal with leaf node
        if(target_node.is_leaf)
        {
            for(const auto& prim_idx: target_node.primitive_index_buffer)
            {
                auto result = intersectEquation( tree.primitive(prim_idx), ray);
                auto hit_t = result(0);
                if(!validIntersect (result) || !ray.valid(hit_t))
                    continue;

                if(hit_type == eAnyHit
                || hit_type == eMultiHit
                || (hit_type == eClosestHit && data_hit_t.empty()))
                {
                    data_hit_t.push_back(hit_t);
                    data_prim_idx.push_back(prim_idx);
                    data_vertex_coeff.push_back(1.f - mxm::sum(result) + hit_t);
                    for(size_t i = 1; i < result.size(); i++)
                        data_vertex_coeff.push_back(result(i));
                }else if (hit_type == eClosestHit && !data_hit_t.empty() && hit_t < data_hit_t.at(0))
                {
                    data_hit_t.at(0) = hit_t;
                    data_prim_idx.at(0) = prim_idx;
                    data_vertex_coeff.at(0) = (1.f - mxm::sum(result) + hit_t);
                    for(size_t i = 1; i < result.size(); i++)
                        data_vertex_coeff.at(i) = result(i);
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

    if(p_coeff) (*p_coeff) = std::move(Matrix<float>({tree.dim(), data_hit_t.size()}, std::move(data_vertex_coeff), COL));
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
    std::stack<Node> stk;
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
                if(distance(nodeBuffer().at(idx).aabb, pt).at(0) < radius)
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
                auto dist = distance(nodeBuffer().at(idx).aabb, pt)[0];
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
