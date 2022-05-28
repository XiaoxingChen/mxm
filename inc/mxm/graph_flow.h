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
    return ret;
}

template<typename WeightType> BinaryEdgeProperty<WeightType, true>
residualCapacityWeight(const BinaryEdgeProperty<WeightType, true>& weight)
{
    BinaryEdgeProperty<WeightType, true> ret;
    size_t edge_num = weight.properties().size();
    Vector<WeightType> weight_buffer(edge_num * 2);

    for(size_t i = 0;i < edge_num; i++)
    {
        weight_buffer(i) = weight.property(i);
        weight_buffer(edge_num + i) = WeightType(0.);
    }
    // return ret.setInvalidProperty(WeightType(0.)).initProperty(weight_buffer);
    ret.setInvalidProperty(WeightType(0.));
    ret.initProperty(weight_buffer);
    return ret;
}

template<typename WeightType> BinaryEdgeProperty<WeightType, true>
flowFromCapacityResidual(
    const BinaryEdgeProperty<WeightType, true>& cap,
    const BinaryEdgeProperty<WeightType, true>& residual)
{
    BinaryEdgeProperty<WeightType, true> ret;
    size_t edge_num = cap.properties().size();
    ret.initProperty(cap.properties() - residual.properties()(Block({0, edge_num}, {})));
    ret.setInvalidProperty(0);
    return ret;
}

template<typename GraphType, typename WeightType>
void noneZeroDfsAnyHit(
    const GraphType& g,
    const BinaryEdgeProperty<WeightType, true>& weight,
    size_t src, size_t dst, std::set<size_t>& visited, std::vector<size_t>& path)
{
    visited.insert(src);
    for(auto & adj: g.adjacency(src))
    {
        if(visited.count(adj) > 0) continue;
        if(weight(src, adj) == 0) continue;
        path.push_back(adj);
        if(adj == dst)
        {
            break;
        }
        noneZeroDfsAnyHit(g, weight, adj, dst, visited, path);
        if(path.back() == dst) break;
        path.pop_back();
    }
    return;
}

template<typename GraphType, typename WeightType>
BinaryEdgeProperty<WeightType, true>
fordFulkersonMaxFLow(
    const GraphType& capacity_graph,
    const BinaryEdgeProperty<WeightType, true>& weight,
    size_t src, size_t sink)
{
    auto operation_g = residualCapacityGraph(capacity_graph);
    auto operation_w = residualCapacityWeight(weight);

    operation_w.setTopology(&operation_g);

    WeightType max_flow = WeightType(0.);
    while(true)
    {
        std::vector<size_t> path{src};
        std::set<size_t> visited;
        noneZeroDfsAnyHit(operation_g, operation_w, src, sink, visited, path);
        if(1 == path.size()) break;
        // std::cout << mxm::to_string(path) << std::endl;
        WeightType bottle_neck = operation_w(path.at(0), path.at(1));
        for(size_t i = 1; i < path.size() - 1; i++)
        {
            bottle_neck = std::min(bottle_neck, operation_w(path.at(i), path.at(i+1)));
            max_flow += bottle_neck;
        }
        for(size_t i = 0; i < path.size() - 1; i++)
        {
            operation_w(path.at(i), path.at(i+1)) -= bottle_neck;
            operation_w(path.at(i+1), path.at(i)) += bottle_neck;
        }
    }
    // std::cout << "max_flow: " << max_flow << std::endl;
    return operation_w;
}

} // namespace mxm


#endif // __GRAPH_FLOW_H__
