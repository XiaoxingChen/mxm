#if !defined(__GRAPH_UTILS_H__)
#define __GRAPH_UTILS_H__

#include "mxm/graph_base.h"
#include <set>

namespace mxm
{

template<typename GraphType>
void connectedComponentsDFS(
    const GraphType& g,
    size_t vertex_idx,
    std::set<size_t>& unvisited,
    std::set<size_t>& s)
{
    if(unvisited.count(vertex_idx) == 0) return ;
    s.insert(vertex_idx);
    unvisited.erase(vertex_idx);
    for(const auto & adj_idx: g.adjacency(vertex_idx))
    {
        connectedComponentsDFS(g, adj_idx, unvisited, s);
    }
}

//
// GraphType requirements:
// size_t GraphType::vertexNum()
// std::vector<size_t> GraphType::adjacency()
template<typename GraphType>
std::vector<std::set<size_t>> connectedComponents(const GraphType& g)
{
    std::set<size_t> unvisited;
    std::vector<std::set<size_t>> ret;
    for(size_t i = 0; i < g.vertexNum(); i++) unvisited.insert(i);

    while(! unvisited.empty())
    {
        ret.push_back(std::set<size_t>());
        auto target = *unvisited.begin();
        connectedComponentsDFS(g, target, unvisited, ret.back());
    }
    return ret;
}

// For undirected graph
template <typename GraphType>
bool isConnected(const GraphType& g)
{
    return connectedComponents(g).size() == 1;
}

template <typename GraphType>
bool isCyclic(const GraphType& g, bool do_connectivity_test=true)
{
    if(do_connectivity_test && !isConnected(g))
        return false;
    return g.vertexNum() + 1 < g.edgeNum();
}

} // namespace mxm


#endif // __GRAPH_UTILS_H__
