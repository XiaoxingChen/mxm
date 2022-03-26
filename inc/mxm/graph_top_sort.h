#if !defined(__GRAPH_TOP_SORT_H__)
#define __GRAPH_TOP_SORT_H__

#include "mxm/graph_base.h"
#include <set>
#include <queue>

namespace mxm
{

// post order
template<typename GraphType>
void topSortDFS(
    const GraphType& g,
    size_t node,
    std::set<size_t>& unvisited,
    std::vector<size_t>& order)
{
    if(unvisited.count(node) == 0) return; //return if visited

    // std::cout << node << std::endl;
    for(const auto & succ : g.successors(node))
    {
        topSortDFS(g, succ, unvisited, order);
    }
    order.push_back(node);
    // std::cout << mxm::to_string(order) << std::endl;

    unvisited.erase(node);
    // std::cout << "finish: " << node << std::endl;
}

template<typename GraphType>
std::vector<size_t> topSort(const GraphType& g)
{
    std::set<size_t> unvisited;
    std::vector<size_t> top_ordering;
    for(size_t i = 0; i < g.vertexNum(); i++) unvisited.insert(i);

    while(!unvisited.empty())
    {
        topSortDFS(g, *unvisited.begin(), unvisited, top_ordering);
    }
    std::reverse(top_ordering.begin(), top_ordering.end());
    return top_ordering;
}

} // namespace mxm



#endif // __GRAPH_TOP_SORT_H__

