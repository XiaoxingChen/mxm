#if !defined(__GRAPH_UTILS_H__)
#define __GRAPH_UTILS_H__

#include "mxm/graph_base.h"
#include <set>

namespace mxm
{

template<typename GraphType>
std::enable_if_t<!GraphType::directed(), void>
connectedComponentsDFS(
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


template<typename GraphType>
std::enable_if_t<!GraphType::directed(), std::vector<std::set<size_t>>>
connectedComponents(const GraphType& g)
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

template<typename GraphType>
void weaklyConnectedComponentsDFS(
    const GraphType& g,
    size_t target,
    std::set<size_t>& unvisited,
    std::vector<std::set<size_t>>& component_list,
    std::set<size_t>& touched_components)
{
    if(unvisited.count(target) == 0)
    {
        // target visited
        for(size_t i = 0; i < component_list.size(); i++)
        {
            if(component_list.at(i).count(target) > 0)
            {
                touched_components.insert(i);
                break;
            }
        }
        return ;
    }
    unvisited.erase(target);
    component_list.back().insert(target);

    for(const auto & adj: g.adjacency(target))
    {
        weaklyConnectedComponentsDFS(g, adj, unvisited, component_list, touched_components);
    }
}

template<typename GraphType>
std::vector<std::set<size_t>>
weaklyConnectedComponents(const GraphType& g)
{
    std::set<size_t> unvisited;
    std::vector<std::set<size_t>> component_list;
    for(size_t i = 0; i < g.vertexNum(); i++) unvisited.insert(i);

    while(! unvisited.empty())
    {
        component_list.push_back(std::set<size_t>());
        auto target = *unvisited.begin();
        std::set<size_t> touched_components; // ascending order
        weaklyConnectedComponentsDFS(g, target, unvisited, component_list, touched_components);

        if(touched_components.size() == 0) continue;
        // merge touched components
        auto min_comp_idx = *std::min_element(touched_components.begin(), touched_components.end());
        touched_components.erase(min_comp_idx);
        for(const auto & comp_idx: touched_components)
        {
            const auto & comp_to_be_merged = component_list.at(comp_idx);
            component_list.at(min_comp_idx).insert(comp_to_be_merged.cbegin(), comp_to_be_merged.cend());
        }
        for(auto it = touched_components.rbegin(); it != touched_components.rend(); it++)
        {
            component_list.erase(component_list.begin() + (*it));
        }
    }
    return component_list;
}

// For undirected graph
template <typename GraphType>
bool isConnected(const GraphType& g)
{
    return connectedComponents(g).size() == 1;
}

template <typename GraphType>
bool isWeaklyConnected(const GraphType& g)
{
    return weaklyConnectedComponents(g).size() == 1;
}

template <typename GraphType>
std::enable_if_t<!GraphType::directed(), bool>
isCyclic(const GraphType& g)
{
    size_t component_num = connectedComponents(g).size();
    return g.vertexNum() < g.edgeNum() + component_num;
}

template <typename GraphType>
std::enable_if_t<GraphType::directed(), bool>
isCyclic(const GraphType& g)
{
    size_t component_num = weaklyConnectedComponents(g).size();
    return g.vertexNum() < g.edgeNum() + component_num;
}

} // namespace mxm


#endif // __GRAPH_UTILS_H__
