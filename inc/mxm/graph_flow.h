#if !defined(__GRAPH_FLOW_H__)
#define __GRAPH_FLOW_H__

#include "graph_base.h"
#include <set>

namespace mxm
{

template<typename GraphType> GraphType
residualCapacityGraph(const GraphType& g)
{
    GraphType ret(g.vertexNum());
    Matrix<size_t> edge_buffer({2, g.edgeNum()*2});
    for(auto & e_idx_pair: g.edgeIndices())
    {
        size_t idx = e_idx_pair.second;
        edge_buffer(0, idx) = e_idx_pair.first[0];
        edge_buffer(1, idx) = e_idx_pair.first[1];
        edge_buffer(0, idx + g.edgeNum()) = e_idx_pair.first[1];
        edge_buffer(1, idx + g.edgeNum()) = e_idx_pair.first[0];
    }
    ret.initEdges(edge_buffer);
    using DType = typename GraphType::WeightType;
    Vector<DType> properties(g.edgeNum()*2);

    for(size_t i = 0;i < g.edgeNum(); i++)
    {
        properties(i) = g.property(i);
        properties(g.edgeNum() + i) = DType(0.);
    }

    ret.initProperty(g.properties());
    return ret;
}

template<typename GraphType>
void noneZeroDfsAnyHit(const GraphType& g, size_t src, size_t dst, std::set<size_t>& visited, std::vector<size_t>& path)
{
    visited.insert(src);
    for(auto & adj: g.adjacency(src))
    {
        if(visited.count(adj) > 0) continue;
        if(g.weight(src, adj) == 0) continue;
        path.push_back(adj);
        if(adj == dst)
        {
            break;
        }
        noneZeroDfsAnyHit(g, adj, dst, visited, path);
        if(path.back() == dst) break;
        path.pop_back();
    }
    return;
}

template<typename GraphType>
GraphType fordFulkersonMaxFLow(const GraphType& capacity_graph, size_t src, size_t sink)
{
    using DType = typename GraphType::WeightType;
    auto operation_g = residualCapacityGraph(capacity_graph);

    DType max_flow = DType(0.);
    while(true)
    {
        std::vector<size_t> path{src};
        std::set<size_t> visited;
        noneZeroDfsAnyHit(operation_g, src, sink, visited, path);
        if(1 == path.size()) break;
        // std::cout << mxm::to_string(path) << std::endl;
        DType bottle_neck = operation_g.weight(path.at(0), path.at(1));
        for(size_t i = 1; i < path.size() - 1; i++)
        {
            bottle_neck = std::min(bottle_neck, operation_g.weight(path.at(i), path.at(i+1)));
            max_flow += bottle_neck;
        }
        for(size_t i = 0; i < path.size() - 1; i++)
        {
            operation_g.property(path.at(i), path.at(i+1)) -= bottle_neck;
            operation_g.property(path.at(i+1), path.at(i)) += bottle_neck;
        }
    }
    // std::cout << "max_flow: " << max_flow << std::endl;
    return operation_g;
}

} // namespace mxm


#endif // __GRAPH_FLOW_H__
