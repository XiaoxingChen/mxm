#if !defined(_SPATIAL_METRIC_TREE_H_)
#define _SPATIAL_METRIC_TREE_H_

#include <vector>
#include <stack>
#include <map>
#include <array>
#include <memory>
#include <chrono>
#include "spatial_tree_base.h"

// References: https://codeforces.com/blog/entry/17974

namespace mxm
{
#if 0
// template<typename ArithType, typename PointType>
// using DistanceFuncType = std::function< ArithType(const PointType&, const PointType&)>;

// template<typename PointType>
// using DistanceFuncType = std::function<typename Traits<PointType>::ArithType (const PointType&, const PointType&)>;
template<typename PointType>
typename Traits<PointType>::ArithType
distance(const PointType&, const PointType&);

// template<typename PointType>
// using MidPointFuncType = std::function< PointType(const PointType&, const PointType&)>;
template<typename PointType>
PointType midPoint(const PointType&, const PointType&);

template<typename PointType>
PointType interpolate(const PointType&, const PointType&, typename Traits<PointType>::ArithType t);

template<typename PointSetType>
typename Traits<PointSetType>::PointType
getPoint(const PointSetType& point_set, size_t i);

template<typename PointSetType>
size_t getSize(const PointSetType& point_set);

template<typename PointSetType>
size_t farthestPoint(
    const PointSetType& point_set,
    size_t target_idx,
    std::vector<size_t>::const_iterator idx_iter_begin,
    std::vector<size_t>::const_iterator idx_iter_end)
{
    auto target_pt = getPoint(point_set, target_idx);
    using ArithType = typename Traits<decltype(target_pt)>::ArithType;
    ArithType max_dist = distance(target_pt, target_pt);
    size_t ret_idx = target_idx;

    for(auto iter = idx_iter_begin; iter != idx_iter_end; iter++)
    {
        auto tmp_pt = getPoint(point_set, *iter);
        ArithType tmp_dist = distance(tmp_pt, target_idx);
        if(tmp_dist > max_dist)
        {
            ret_idx = *iter;
            max_dist = tmp_dist;
        }
    }
    return ret_idx;
}

template<typename PointSetType> std::array<size_t, 2>
farthestPair(
    const PointSetType& point_set,
    std::vector<size_t>::const_iterator idx_iter_begin,
    std::vector<size_t>::const_iterator idx_iter_end)
{
    size_t p0 = farthestPoint(point_set, *idx_iter_begin, idx_iter_begin, idx_iter_end);
    size_t p1 = farthestPoint(point_set, p0, idx_iter_begin, idx_iter_end);
    return {p0, p1};
}

template<typename PointType>
bool biasRate(const PointType& target, const PointType& p0, const PointType& p1)
{
    // the closer to p1, the larger
    auto dist_p0 = distance(target, p0);
    auto dist_p1 = distance(target, p1);
    return dist_p0 - dist_p1;
}

template<typename PointSetType>
class BallTree
{
public:
    using ThisType = BallTree<PointSetType>;
    using PointType = typename Traits<PointSetType>::PointType;
    using ArithType = typename Traits<PointType>::ArithType;
    constexpr size_t pointPerLeaf() const {return 2;}
    struct Node
    {
        Node():is_leaf(false){}
        PointType center;
        ArithType radius;
        std::vector<size_t> children_index_buffer;
        std::vector<size_t> point_index_buffer;
        bool is_leaf;
        size_t point_num;
    };
    void build(bool verbose=true);
private:
    void updateBoundary(size_t root_node_idx);
    std::vector<ThisType::Node> node_buffer_;
    std::shared_ptr<PointSetType> point_set_;
};

size_t treeNodeRequirement(size_t);

template<typename PointSetType>
void BallTree<PointSetType>::updateBoundary(size_t root_node_idx)
{
    ThisType::Node& target_node = node_buffer_.at(root_node_idx);
    if(target_node.is_leaf)
    {
        if(1 == target_node.point_num)
        {
            target_node.center = getPoint(point_set_, target_node.point_index_buffer.at(0));
            target_node.radius = 0.;
        }else if(2 == target_node.point_num)
        {
            auto p0 = getPoint(point_set_, target_node.point_index_buffer.at(0));
            auto p1 = getPoint(point_set_, target_node.point_index_buffer.at(1));
            target_node.center = interpolate(p0, p1, 0.5);
            target_node.radius = distance(p0, p1);
        }
    }else //internal node
    {
        for(size_t i = 0; i < target_node.children_index_buffer.size(); i++)
        {
            auto child_idx = target_node.children_index_buffer.at(i);
            updateBoundary(child_idx);
        }
        auto child0 = node_buffer_.at(target_node.children_index_buffer.at(0));
        auto child1 = node_buffer_.at(target_node.children_index_buffer.at(1));
        target_node.radius = (distance(child0.center, child1.center) + child0.radius + child1.radius) * 0.5;
        // target_node.center = interpolate(child0.center, child1.center, 0.5); //todo
    }
}
#endif

size_t treeNodeRequirement(size_t leaf_num)
{
    const size_t child_num = 2;
    size_t height = static_cast<size_t>(1. + std::log(leaf_num) / std::log(child_num)) + 1;
    return (1 << height);
}

template<typename PointSetType, typename PointType, typename ArithType>
class MetricTreeBase
{
public:

    struct Sphere
    {
        size_t center_idx;
        ArithType radius;
    };

    struct Node
    {
        Node():is_leaf(false){}
        Sphere bounding;
        std::vector<size_t> children_index_buffer;
        std::vector<size_t> point_index_buffer;
        bool is_leaf;
        size_t point_num;
    };

    using ThisType = MetricTreeBase<PointSetType, PointType, ArithType>;

    virtual PointType point(size_t idx) const = 0;
    virtual ArithType distance(const PointType& p0, const PointType& p1) const = 0;
    virtual size_t size() const = 0;
    void build(bool verbose=true);
    constexpr size_t pointPerLeaf() const {return 2;}

    std::multimap<ArithType, size_t> radiusSearch(const PointType& p0, ArithType radius) const;
    std::multimap<ArithType, size_t> nearestNeighborSearch(const PointType& p0, size_t k) const;


    bool verifyTree(size_t node_idx)
    {
        if(node_buffer_.at(node_idx).is_leaf)
        {
            for(auto pt_idx: node_buffer_.at(node_idx).point_index_buffer)
            {
                if(fastDistance(pt_idx, node_buffer_.at(node_idx).bounding.center_idx) > node_buffer_.at(node_idx).bounding.radius)
                {
                    return false;
                }
            }
            return true;
        }

        for(auto & child_idx: node_buffer_.at(node_idx).children_index_buffer)
        {
            if(fastDistance(
                node_buffer_.at(child_idx).bounding.center_idx,
                node_buffer_.at(node_idx).bounding.center_idx)
                > node_buffer_.at(node_idx).bounding.radius)
                return false;
            if(!verifyTree(child_idx))
                return false;
        }
        return true;
    }

protected:
    std::map<std::array<size_t, 2>, ArithType> distCache_;
    std::vector<ThisType::Node> node_buffer_;



protected:
    std::map<size_t, ArithType> updateBoundary(size_t root_node_idx);
    ArithType biasRate(size_t target, size_t p0, size_t p1)
    {
        // the closer to p1, the larger
        auto dist_p0 = fastDistance(target, p0);
        auto dist_p1 = fastDistance(target, p1);
        return dist_p0 - dist_p1;
    }

    ArithType fastDistance(size_t p0, size_t p1)
    {
        if(p0 == p1) return ArithType(0.);
        std::array<size_t, 2> key = p0 < p1 ? std::array<size_t, 2>{p0, p1} : std::array<size_t, 2>{p1, p0};
        if(distCache_.count(key) == 0)
            distCache_[key] = distance(point(p0), point(p1));
        return distCache_[key];
    }

    // complexity: ?
    std::map<size_t, ArithType> mergeEccentricity(
        const std::map<size_t, ArithType>& eccent_0,
        const std::map<size_t, ArithType>& eccent_1)
    {
        std::map<size_t, ArithType> ret_eccent;
        for(const auto& pair : eccent_0) ret_eccent.insert(pair);
        for(const auto& pair : eccent_1) ret_eccent.insert(pair);

        for(const auto& idx_radius_pair0 : eccent_0)
        {
            size_t idx0 = idx_radius_pair0.first;
            for(const auto& idx_radius_pair1 : eccent_1)
            {
                size_t idx1 = idx_radius_pair1.first;
                auto dist = fastDistance(idx0, idx1);
                ret_eccent[idx0] = std::max(ret_eccent[idx0], dist);
                ret_eccent[idx1] = std::max(ret_eccent[idx1], dist);
            }
        }
        return ret_eccent;
    }

#if 1
    // complexity: O(n)
    ArithType eccentricity(
        size_t target_idx,
        std::vector<size_t>::const_iterator idx_iter_begin,
        std::vector<size_t>::const_iterator idx_iter_end)
    {
        auto ret_idx = farthestPoint(target_idx, idx_iter_begin, idx_iter_end);
        return fastDistance(target_idx, ret_idx);
    }

    // complexity: O(n^2)
    std::map<size_t, ArithType> rangeEccentricity(
        std::vector<size_t>::const_iterator idx_iter_begin,
        std::vector<size_t>::const_iterator idx_iter_end)
    {
        std::map<size_t, ArithType> ret;
        for(auto iter = idx_iter_begin; iter != idx_iter_end; iter++)
        {
            ret[*iter] = eccentricity(*iter, idx_iter_begin, idx_iter_end);
        }
        return ret;
    }
#endif
    // complexity: O(n)
    size_t farthestPoint(
        size_t target_idx,
        std::vector<size_t>::const_iterator idx_iter_begin,
        std::vector<size_t>::const_iterator idx_iter_end)
    {
        ArithType max_dist(0);
        size_t ret_idx = target_idx;

        for(auto iter = idx_iter_begin; iter != idx_iter_end; iter++)
        {
            ArithType tmp_dist = fastDistance(*iter, target_idx);
            if(tmp_dist > max_dist)
            {
                ret_idx = *iter;
                max_dist = tmp_dist;
            }
        }
        return ret_idx;
    }

    std::array<size_t, 2>
    farthestPair(
        std::vector<size_t>::const_iterator idx_iter_begin,
        std::vector<size_t>::const_iterator idx_iter_end)
    {
        size_t p0 = farthestPoint(*idx_iter_begin, idx_iter_begin, idx_iter_end);
        size_t p1 = farthestPoint(p0, idx_iter_begin, idx_iter_end);
        return {p0, p1};
    }
};

template<typename PointSetType, typename PointType, typename ArithType>
void MetricTreeBase<PointSetType, PointType, ArithType>::build(bool verbose)
{
    std::vector<size_t> point_index_buffer(size());
    for(size_t i = 0; i < point_index_buffer.size(); i++) point_index_buffer.at(i) = i;

    {
        size_t require = treeNodeRequirement(point_index_buffer.size() / pointPerLeaf() + 1);
        node_buffer_.clear();
        node_buffer_.reserve(require);
    }

    // DFS Preorder
    std::stack<RangeNode> stk;

    stk.push({node_buffer_.size(), {0, point_index_buffer.size()}});
    node_buffer_.push_back( ThisType::Node() );

    auto t_start = std::chrono::system_clock::now();

    while (!stk.empty())
    {
        auto target = stk.top();
        stk.pop();

        auto& target_node = node_buffer_.at(target.node_idx);
        target_node.point_num = target.range[1] - target.range[0];
        if(target_node.point_num <= pointPerLeaf())
        {
            target_node.is_leaf = true;
            for(size_t i = target.range[0]; i < target.range[1]; i++)
            {
                target_node.point_index_buffer.push_back(point_index_buffer.at(i));
            }
        }else // is internal node
        {
            auto polar_pair = farthestPair(
                point_index_buffer.begin() + target.range[0],
                point_index_buffer.begin() + target.range[1]);

            size_t mid = (target.range[1] + target.range[0]) / 2;

            std::nth_element(
                point_index_buffer.begin() + target.range[0],
                point_index_buffer.begin() + mid,
                point_index_buffer.begin() + target.range[1],
                [&](const size_t& point_idx1, const size_t& point_idx2)
                {
                    return biasRate(point_idx1, polar_pair[0], polar_pair[1]) < biasRate(point_idx2, polar_pair[0], polar_pair[1]);
                }
            );

            // create children
            stk.push({node_buffer_.size(), {target.range[0], mid}});
            target_node.children_index_buffer.push_back(node_buffer_.size());
            node_buffer_.push_back( ThisType::Node() );

            stk.push({node_buffer_.size(), {mid, target.range[1]}});
            target_node.children_index_buffer.push_back(node_buffer_.size());
            node_buffer_.push_back( ThisType::Node() );
        }
    }
    auto t_end = std::chrono::system_clock::now();
    std::cout << "build tree t_cost(s): " << std::chrono::duration<double>(t_end - t_start).count() << std::endl;

    updateBoundary(0);
}

template<typename PointSetType, typename PointType, typename ArithType>
std::map<size_t, ArithType>
MetricTreeBase<PointSetType, PointType, ArithType>::updateBoundary(size_t root_node_idx)
{
    ThisType::Node& target_node = node_buffer_.at(root_node_idx);
    std::map<size_t, ArithType> ret_eccentricity;
    if(target_node.is_leaf)
    {
        ret_eccentricity = rangeEccentricity(target_node.point_index_buffer.begin(), target_node.point_index_buffer.end());

    }else //internal node
    {
        auto local_eccent_0 = updateBoundary(target_node.children_index_buffer.at(0));
        auto local_eccent_1 = updateBoundary(target_node.children_index_buffer.at(1));

        ret_eccentricity = mergeEccentricity(local_eccent_0, local_eccent_1);
    }
    auto it = std::min_element(ret_eccentricity.begin(), ret_eccentricity.end(),
            [](const auto& p0, const auto& p1){return p0.second < p1.second;});
    target_node.bounding.center_idx = it->first;
    target_node.bounding.radius = it->second;

    return ret_eccentricity;
}

template<typename PointSetType, typename PointType, typename ArithType>
std::multimap<ArithType, size_t>
MetricTreeBase<PointSetType, PointType, ArithType>::nearestNeighborSearch(const PointType& p0, size_t k) const
{
    std::multimap<float, size_t> ret;

    std::stack<size_t> stk;
    stk.push(0);
    while(!stk.empty())
    {
        auto target = node_buffer_.at(stk.top());
        stk.pop();

        auto curr_max_dist = ret.empty() ? INFINITY : std::prev(ret.end())->first;

        if(target.is_leaf)
        {
            for(auto idx: target.point_index_buffer)
            {
                auto dist = distance(point(idx), p0);
                if(dist < curr_max_dist)
                    ret.insert({dist, idx});
            }
        }else
        {
            for(auto node_idx: target.children_index_buffer)
            {
                auto pt_center = point(node_buffer_.at(node_idx).bounding.center_idx);
                auto dist = distance(p0, pt_center);
                if(dist - node_buffer_.at(node_idx).bounding.radius < curr_max_dist)
                {
                    stk.push(node_idx);
                }
            }
        }

        while(ret.size() > k)
        {
            ret.erase(std::prev(ret.end()));
        }

    }

    return ret;
}

template<typename PointSetType, typename PointType, typename ArithType>
std::multimap<ArithType, size_t>
MetricTreeBase<PointSetType, PointType, ArithType>::radiusSearch(const PointType& p0, ArithType radius) const
{
    std::multimap<ArithType, size_t> ret;
    std::stack<size_t> stk;
    stk.push(0);

    while(!stk.empty())
    {
        const ThisType::Node& target_node = node_buffer_.at(stk.top());
        stk.pop();

        if(target_node.is_leaf)
        {
            for(size_t idx: target_node.point_index_buffer)
            {
                auto dist = distance(p0, point(idx));
                if(dist <= radius)
                {
                    ret.insert({dist, idx});
                }
            }
        }else // internal node
        {
            for(size_t node_idx: target_node.children_index_buffer)
            {
                auto pt_center = point(node_buffer_.at(node_idx).bounding.center_idx);
                auto dist = distance(p0, pt_center);
                if(dist - node_buffer_.at(node_idx).bounding.radius < radius)
                {
                    stk.push(node_idx);
                }
            }
        }
    }
    return ret;
}

} // namespace mxm


#endif // _SPATIAL_METRIC_TREE_H_