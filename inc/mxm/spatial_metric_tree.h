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
    };

    using ThisType = MetricTreeBase<PointSetType, PointType, ArithType>;

    virtual PointType point(size_t idx) const = 0;
    virtual ArithType distance(const PointType& p0, const PointType& p1) const = 0;
    virtual size_t size() const = 0;
    std::multimap<ArithType, size_t> build1(size_t root_node_idx, int verbose=0);
    void build(int verbose=0);
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
    std::vector<size_t> index_storage_buffer_;

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
std::string to_string(const std::vector<typename MetricTreeBase<PointSetType, PointType, ArithType>::Node>& node_buffer);


template<typename PointSetType, typename PointType, typename ArithType>
std::multimap<ArithType, size_t>
MetricTreeBase<PointSetType, PointType, ArithType>::build1(size_t root_node_idx, int verbose)
{
    if(0 == root_node_idx)
    {
        index_storage_buffer_.resize(size());
        for(size_t i = 0; i < index_storage_buffer_.size(); i++) index_storage_buffer_.at(i) = i;

        size_t require = treeNodeRequirement(index_storage_buffer_.size() / pointPerLeaf() + 1);
        node_buffer_.clear();
        node_buffer_.reserve(require);
        auto t_start = std::chrono::system_clock::now();
    }

    node_buffer_.push_back( ThisType::Node() );

    auto& target_node = node_buffer_.back();
    return {};

}

template<typename PointSetType, typename PointType, typename ArithType>
void MetricTreeBase<PointSetType, PointType, ArithType>::build(int verbose)
{
    struct Context
    {
        size_t node_idx = 0;
        std::array<size_t, 2> range{0,0};
        size_t child_idx = 0;
        size_t mid_idx = 0;
        std::vector<std::map<size_t, ArithType>> child_eccentricities;
    };

    std::vector<size_t> point_index_buffer(size());
    for(size_t i = 0; i < point_index_buffer.size(); i++) point_index_buffer.at(i) = i;

    {
        size_t require = treeNodeRequirement(point_index_buffer.size() / pointPerLeaf() + 1);
        node_buffer_.clear();
        node_buffer_.reserve(require);
    }

    // DFS
    std::stack<Context> stk;

    stk.push({node_buffer_.size(), {0, point_index_buffer.size()}});
    node_buffer_.push_back( ThisType::Node() );

    auto t_start = std::chrono::system_clock::now();

    while (!stk.empty())
    {
        auto& context = stk.top();

        auto& target_node = node_buffer_.at(context.node_idx);
        auto point_num = context.range[1] - context.range[0];
        if(point_num <= pointPerLeaf())
        {
            target_node.is_leaf = true;
            for(size_t i = context.range[0]; i < context.range[1]; i++)
            {
                target_node.point_index_buffer.push_back(point_index_buffer.at(i));
            }

            auto local_eccentricity = rangeEccentricity(target_node.point_index_buffer.begin(), target_node.point_index_buffer.end());
            auto it = std::min_element(local_eccentricity.begin(), local_eccentricity.end(),
                    [](const auto& p0, const auto& p1){return p0.second < p1.second;});
            target_node.bounding.center_idx = it->first;
            target_node.bounding.radius = it->second;

            stk.pop();
            if(!stk.empty()) // return result to stake top
                stk.top().child_eccentricities.push_back(std::move(local_eccentricity));

            continue; // make stack popping work
        }
        // is internal node

        if(0 == context.child_idx)
        {
            auto polar_pair = farthestPair(
                point_index_buffer.begin() + context.range[0],
                point_index_buffer.begin() + context.range[1]);

            context.mid_idx = (context.range[1] + context.range[0]) / 2;

            std::nth_element(
                point_index_buffer.begin() + context.range[0],
                point_index_buffer.begin() + context.mid_idx,
                point_index_buffer.begin() + context.range[1],
                [&](const size_t& point_idx1, const size_t& point_idx2)
                {
                    return biasRate(point_idx1, polar_pair[0], polar_pair[1]) < biasRate(point_idx2, polar_pair[0], polar_pair[1]);
                }
            );

            // create children
            stk.push({node_buffer_.size(), {context.range[0], context.mid_idx}});
            target_node.children_index_buffer.push_back(node_buffer_.size());
            node_buffer_.push_back( ThisType::Node() );

            context.child_idx++;
            continue; // make stack pushing work
        }else if (1 == context.child_idx)
        {
            stk.push({node_buffer_.size(), {context.mid_idx, context.range[1]}});
            target_node.children_index_buffer.push_back(node_buffer_.size());
            node_buffer_.push_back( ThisType::Node() );

            context.child_idx++;
            continue; // make stack pushing work
        }else // post process
        {
            auto local_eccentricity = mergeEccentricity(
                context.child_eccentricities.at(0),
                context.child_eccentricities.at(1));
            auto it = std::min_element(local_eccentricity.begin(), local_eccentricity.end(),
                [](const auto& p0, const auto& p1){return p0.second < p1.second;});
            target_node.bounding.center_idx = it->first;
            target_node.bounding.radius = it->second;

            stk.pop();
            if(!stk.empty()) // return result to stake top
                stk.top().child_eccentricities.push_back(std::move(local_eccentricity));
        }
    }
    auto t_end = std::chrono::system_clock::now();
    if(verbose > 0)
    {
        std::cout << "build tree t_cost(s): " << std::chrono::duration<double>(t_end - t_start).count() << std::endl;
    }
    if(verbose > 1)
    {
        std::cout << to_string<PointSetType, PointType, ArithType>(node_buffer_) << std::endl;
    }



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

template<typename PointSetType, typename PointType, typename ArithType>
std::string to_string(const std::vector<typename MetricTreeBase<PointSetType, PointType, ArithType>::Node>& node_buffer)
{
    std::string ret;
    for(size_t i = 0; i < node_buffer.size(); i++)
    {
        ret += "{idx: ";
        ret += (std::to_string(i) + ", center_idx: " + std::to_string(node_buffer.at(i).bounding.center_idx));
        ret += (", radius: " + std::to_string(node_buffer.at(i).bounding.radius));

        ret += ", children_index_buffer: [";
        for(size_t j=0; j < node_buffer.at(i).children_index_buffer.size(); j++)
        {
            size_t child_idx = node_buffer.at(i).children_index_buffer.at(j);
            ret += std::to_string(child_idx);
            if(j != node_buffer.at(i).children_index_buffer.size()-1) ret += ",";
        } ret += "]";

        ret += ", point_index_buffer: [";
        for(size_t j=0; j < node_buffer.at(i).point_index_buffer.size(); j++)
        {
            size_t child_idx = node_buffer.at(i).point_index_buffer.at(j);
            ret += std::to_string(child_idx);
            if(j != node_buffer.at(i).point_index_buffer.size()-1) ret += ",";
        } ret += "]";
        ret += "}";
        if(i != node_buffer.size() - 1) ret += ",";

    }

    return ret;
}

} // namespace mxm


#endif // _SPATIAL_METRIC_TREE_H_