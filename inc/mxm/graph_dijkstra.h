#if !defined(_GRAPH_DIJKSTRA_H)
#define _GRAPH_DIJKSTRA_H

#include "graph_base.h"
#include <set>

namespace mxm
{

template<typename GraphType>
std::vector<size_t> dijkstra(
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
            DType local_distance = g.distance(target_idx, succ_idx);

            if(local_distance < g.distance(best_pred.at(succ_idx), succ_idx))
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
typename GraphType::DistanceType
pathFromBestPredecessor(
    const GraphType& g,
    const std::vector<size_t>& best_pred,
    size_t destination,
    std::vector<size_t>* p_path=nullptr)
{
    if(p_path) p_path->clear();
    typename GraphType::DistanceType distance = 0;

    size_t p = destination;
    while(best_pred.at(p) != p)
    {
        size_t pred = best_pred.at(p);
        if(! g.validVertex(pred))
        {
            if(p_path) p_path->clear();
            return INFINITY;
        }
        distance += g.distance(pred, p);
        p = pred;
        if(p_path) p_path->push_back(p);
    }

    if(p_path)
    {
        std::reverse(p_path->begin(), p_path->end());
        p_path->push_back(destination);
    }

    return distance;
}

} // namespace mxm


#endif // _GRAPH_DIJKSTRA_H
