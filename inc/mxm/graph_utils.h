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
    for(const auto & neighbor_idx: g.neighbors(vertex_idx))
    {
        connectedComponentsDFS(g, neighbor_idx, unvisited, s);
    }
}

//
// GraphType requirements:
// size_t GraphType::vertexNum()
// std::vector<size_t> GraphType::neighbors()
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

} // namespace mxm


#endif // __GRAPH_UTILS_H__
