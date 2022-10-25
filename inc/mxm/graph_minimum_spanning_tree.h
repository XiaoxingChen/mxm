#if !defined(__GRAPH_MINIMUM_SPANNING_TREE_H__)
#define __GRAPH_MINIMUM_SPANNING_TREE_H__

#include "graph_utils.h"

namespace mxm
{


// Reference: 
// https://www.youtube.com/watch?v=JZBQLXgSGfs&list=PLDV1Zeh2NRsBI1C-mR6ZhHTyfoEJWlxvq&index=2
template <typename GraphType, typename WeightType>
std::enable_if_t<!GraphType::directed(), GraphType>
kruskalMinimumSpanningTree(const GraphType& g,
    const BinaryEdgeProperty<WeightType, GraphType::directed()> f_dist)
{
    // create sorted edge
    std::vector<std::pair<WeightType, std::array<size_t, 2> > > sorted_edges;
    sorted_edges.reserve(g.edgeNum());
    for(auto & e_idx_pair: g.edgeIndices())
    {
        sorted_edges.push_back({
            f_dist(e_idx_pair.first[0], e_idx_pair.first[1]),
            e_idx_pair.first});
    }
    std::sort(sorted_edges.begin(), sorted_edges.end());

    // create edge buffer
    size_t edge_cnt = 0;
    Matrix<size_t> mst_edge_buff({2, g.vertexNum() - 1});

    // do union find
    std::vector<size_t> union_find_array(g.edgeNum());
    for(size_t i = 0; i < g.edgeNum(); i++) union_find_array.at(i) = i;

    for(const auto & weight_edge: sorted_edges)
    {
        size_t idx_0 = weight_edge.second[0];
        size_t idx_1 = weight_edge.second[1];
        size_t root_0 = idx_0;
        size_t root_1 = idx_1;
        while(union_find_array.at(root_0) != root_0) root_0 = union_find_array.at(root_0);
        while(union_find_array.at(root_1) != root_1) root_1 = union_find_array.at(root_1);

        if(root_0 == root_1) continue;
        if(root_0 == idx_0)
        {
            union_find_array.at(idx_0) = root_1;
        }else if(root_1 == idx_1)
        {
            union_find_array.at(idx_1) = root_0;
        } else
        {
            union_find_array.at(root_0) = root_1;
        }
        mst_edge_buff(0, edge_cnt) = idx_0;
        mst_edge_buff(1, edge_cnt) = idx_1;
        edge_cnt ++;
    }

    GraphType union_graph(g.vertexNum());
    union_graph.initEdges(mst_edge_buff);
    return union_graph;
}
} // namespace mxm

#endif // __GRAPH_MINIMUM_SPANNING_TREE_H__
