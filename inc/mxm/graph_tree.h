#if !defined(__GRAPH_TREE_H__)
#define __GRAPH_TREE_H__

#include "mxm/graph_base.h"

namespace mxm
{

template <typename GraphType>
void eulerTourDFS(
    const GraphType& g,
    size_t target,
    std::vector<std::array<size_t, 2>>& path,
    size_t curr_depth)
{
    path.push_back({target, curr_depth});
    for(const auto& child: g.adjacency(target))
    {
        eulerTourDFS(g, child, path, curr_depth + 1);
    }
    path.push_back({target, curr_depth});
}

template <typename GraphType>
std::enable_if_t<GraphType::directed(), std::vector<std::array<size_t, 2>> >
eulerTour(const GraphType& g, size_t root)
{
    std::array< std::vector<size_t>, 2> path_with_depth;
    eulerTourDFS(g, root, path_with_depth, 0);
    return path_with_depth;
}

} // namespace mxm



#endif // __GRAPH_TREE_H__
