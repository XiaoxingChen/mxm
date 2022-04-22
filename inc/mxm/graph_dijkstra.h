#if !defined(_GRAPH_DIJKSTRA_H)
#define _GRAPH_DIJKSTRA_H

#include "graph_base.h"
#include <set>
#include <queue>

namespace mxm
{

// Summary:Move vertices from Unvisited Set to Visited Set.
// check every vertex is a better predecessor in the Unvisited Set.
template<typename GraphType>
std::enable_if_t< is_weighted_binary_edge<GraphType>::value , std::vector<size_t>>
dijkstraBestPredecessor(
    const GraphType& g,
    size_t start)
{
    using DType = typename GraphType::DistanceType;
    std::multimap<DType, size_t> candidates;
    candidates.insert({DType(0), start});
    // std::set<size_t> candidates;
    // candidates.insert(start);
    std::vector<DType> dist_to_start(g.vertexNum(), std::numeric_limits<DType>::max());
    // std::unordered_map<size_t, DType> dist_to_start;
    dist_to_start.at(start) = DType(0);


    std::set<size_t> visited;

    std::vector<size_t> best_pred(g.vertexNum(), start);

    while(! candidates.empty())
    {
        size_t target_idx = candidates.begin()->second;
        candidates.erase(candidates.begin());
        if(visited.count(target_idx) > 0) continue;

        for(auto & succ_idx: g.adjacency(target_idx))
        {
            if(visited.count(succ_idx) > 0) continue;
            DType latest_distance = g.weight(target_idx, succ_idx) + dist_to_start.at(target_idx);

            if(latest_distance < dist_to_start.at(succ_idx))
            {
                best_pred.at(succ_idx) = target_idx;
                dist_to_start[succ_idx] = latest_distance;
                candidates.insert({latest_distance, succ_idx});
            }
        }

        // move insert before for loop?
        // Seems not necessary because not using DFS
        visited.insert(target_idx); // todo

    }

    // boundary condition
    best_pred.at(start) = start;
    // std::cout << "best_pred: " << mxm::to_string(best_pred) << std::endl;
    return best_pred;
}

template<typename GraphType>
std::enable_if_t< is_unweighted_binary_edge<GraphType>::value , std::vector<size_t>>
dijkstraBestPredecessor(
    const GraphType& g,
    size_t start)
{
    std::queue<size_t> que;
    que.push(start);
    std::set<size_t> visited;

    std::vector<size_t> best_pred(g.vertexNum(), start);

    while(!que.empty())
    {
        size_t target = que.front();
        que.pop();

        visited.insert(target);

        for(const auto & succ : g.adjacency(target))
        {
            if(visited.count(succ) > 0) continue;
            best_pred.at(succ) = target;
            que.push(succ);
        }
    }

    best_pred.at(start) = start;
    return best_pred;
}

template<typename GraphType>
std::enable_if_t< is_weighted_binary_edge<GraphType>::value , std::vector<size_t>>
pathFromBestPredecessor(
    const GraphType& g,
    const std::vector<size_t>& best_pred,
    size_t destination,
    typename GraphType::DistanceType& distance)
{
    std::vector<size_t> path;
    path.clear();
    distance = 0;

    size_t p = destination;
    while(best_pred.at(p) != p)
    {
        size_t pred = best_pred.at(p);
        if(! g.validVertex(pred))
        {
            path.clear();
            distance = INFINITY;
            return path;
        }
        distance += g.weight(pred, p);
        p = pred;
        path.push_back(p);
    }

    {
        std::reverse(path.begin(), path.end());
        path.push_back(destination);
    }

    return path;
}

template<typename GraphType>
std::enable_if_t< is_unweighted_binary_edge<GraphType>::value , std::vector<size_t>>
pathFromBestPredecessor(
    const GraphType& g,
    const std::vector<size_t>& best_pred,
    size_t destination,
    int& distance)
{
    std::vector<size_t> path;
    path.push_back(destination);
    distance = 0;

    auto iterator = destination;
    while(best_pred.at(iterator) != iterator)
    {
        if(!g.validVertex(best_pred.at(iterator)))
        {
            distance = -1;
            return {};
        }
        path.push_back(best_pred.at(iterator));
        distance++;
    }

    return path;
}

template<typename GraphType>
std::enable_if_t<is_weighted_binary_edge<GraphType>::value , std::vector<size_t>>
dijkstra(
    const GraphType& g,
    size_t start,
    size_t dest,
    typename GraphType::DistanceType* p_distance)
{
    auto best_predecessor = dijkstraBestPredecessor(g, start);
    // std::cout << "best_predecessor: " << mxm::to_string(best_predecessor) << std::endl;
    typename GraphType::DistanceType distance;
    auto path = pathFromBestPredecessor(g, best_predecessor, dest, distance);
    if(p_distance) *p_distance = distance;
    return path;
}

template<typename GraphType>
std::enable_if_t< is_unweighted_binary_edge<GraphType>::value , std::vector<size_t>>
dijkstra(
    const GraphType& g,
    size_t start,
    size_t dest,
    int* p_distance)
{
    auto best_predecessor = dijkstraBestPredecessor(g, start);
    typename GraphType::DistanceType distance;
    auto path = pathFromBestPredecessor(g, best_predecessor, dest, distance);
    if(p_distance) *p_distance = distance;
    return path;
}

} // namespace mxm


#endif // _GRAPH_DIJKSTRA_H
