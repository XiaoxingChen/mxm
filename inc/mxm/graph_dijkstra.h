#if !defined(_GRAPH_DIJKSTRA_H)
#define _GRAPH_DIJKSTRA_H

#include "graph_base.h"
#include <set>

namespace mxm
{

inline std::vector<size_t> dijkstra(const SparseDirectedWeightedGraph& g, size_t start)
{
    std::multimap<FloatType, size_t> que;
    que.insert({0, start});

    std::set<size_t> visited;
    visited.insert(start);

    std::vector<size_t> best_pred(g.nodeNum(), g.nodeNum());

    while(! que.empty())
    {
        size_t target_idx = que.begin()->second;
        que.erase(que.begin());

        for(auto & dist_idx_pair: g.successors(target_idx))
        {
            size_t succ_idx = dist_idx_pair.second;
            if(visited.count(succ_idx) > 0) continue;

            if(g.distance(target_idx, succ_idx) < g.distance(best_pred.at(succ_idx), succ_idx))
            {
                best_pred.at(succ_idx) = target_idx;
            }

            que.insert(dist_idx_pair);
        }
        visited.insert(target_idx);
    }

    best_pred.at(start) = start;
    return best_pred;
}

FloatType pathFromBestPredecessor(
    const SparseDirectedWeightedGraph& g,
    const std::vector<size_t>& best_pred,
    size_t destination,
    std::vector<size_t>& path)
{
    path.clear();
    FloatType distance = 0;

    size_t p = destination;
    while(best_pred.at(p) != p)
    {
        size_t pred = best_pred.at(p);
        if(! g.valid(pred))
        {
            path.clear();
            return INFINITY;
        }
        distance += g.distance(pred, p);
        p = pred;
        path.push_back(p);
    }
    std::reverse(path.begin(), path.end());
    path.push_back(destination);
    return distance;
}

} // namespace mxm


#endif // _GRAPH_DIJKSTRA_H
