#if !defined(_GRAPH_DIJKSTRA_H)
#define _GRAPH_DIJKSTRA_H

#include "graph_base.h"
#include <set>

namespace mxm
{

template<typename GraphType>
std::vector<size_t> dijkstraBestPredecessor(
    const GraphType& g,
    size_t start)
{
    using DType = typename GraphType::DistanceType;
    std::multimap<DType, size_t> que;
    que.insert({0, start});

    std::set<size_t> visited;
    visited.insert(start);

    std::vector<size_t> best_pred(g.vertexNum(), g.vertexNum());

    while(! que.empty())
    {
        size_t target_idx = que.begin()->second;
        que.erase(que.begin());

        for(auto & succ_idx: g.successors(target_idx))
        {
            if(visited.count(succ_idx) > 0) continue;
            DType local_distance = g.weight(target_idx, succ_idx);

            if(local_distance < g.weight(best_pred.at(succ_idx), succ_idx))
            {
                best_pred.at(succ_idx) = target_idx;
            }

            que.insert({local_distance, succ_idx});
        }
        visited.insert(target_idx);
    }

    best_pred.at(start) = start;
    return best_pred;
}

template<typename GraphType>
std::vector<size_t>
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
std::vector<size_t> dijkstra(
    const GraphType& g,
    size_t start,
    size_t dest,
    typename GraphType::DistanceType* p_distance)
{
    auto best_predecessor = dijkstraBestPredecessor(g, start);
    typename GraphType::DistanceType distance;
    auto path = pathFromBestPredecessor(g, best_predecessor, dest, distance);
    if(p_distance) *p_distance = distance;
    return path;
}

} // namespace mxm


#endif // _GRAPH_DIJKSTRA_H
