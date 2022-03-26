#if !defined(_TEST_GRAPH_H_)
#define _TEST_GRAPH_H_

#include "test_config.h"

#if TEST_AVAILABLE_ALL
#include <iostream>
#include "mxm/graph_dijkstra.h"
#include "mxm/graph_top_sort.h"
#endif

using namespace mxm;
#if TEST_AVAILABLE_ALL
// test graph:
// link: https://en.wikipedia.org/wiki/File:Dijkstra_Animation.gif
inline void testDijkstra1()
{
    WeightedDirectedGraph<float> g(7);
    Matrix<size_t> edge_buffer(fixRow(2), {1,2, 1,3, 1,6, 2,4, 2,3, 3,6, 3,4, 6,5, 4,5}, COL);
    Vector<float> weight_buffer{7,9,14,15,10,2,11,9,6};
    g.initEdges(edge_buffer);
    g.initWeight(weight_buffer);

    auto best_pred = dijkstra(g, 1);
    std::vector<size_t> path;
    auto distance = pathFromBestPredecessor(g, best_pred, 5, &path);

    if(path != std::vector<size_t>{1,3,6,5})
    {
        std::cout << "best_pred: " << mxm::to_string(best_pred) << std::endl;
        std::cout << "path: " << mxm::to_string(path) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(distance != 20)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

template<typename GraphType>
bool validateTopologicalOrder(const GraphType& g, const std::vector<size_t>& order )
{
    std::set<size_t> built;
    for(size_t i = 0; i < order.size(); i++)
    {
        auto target = order.at(i);
        for(const auto& succ: g.successors(target))
        {
            if(built.count(succ) > 0)
            {
                std::cout << "target: " << target << ", succ: " << succ << std::endl;
                return false;
            }
        }
    }
    return true;
}

inline void testTopologicalSort()
{
    WeightedDirectedGraph<float> g(13);
    Matrix<size_t> edge_buffer(fixRow(2), {
        0,3, 1,3, 2,0, 2,1, 3,6, 3,7, 4,0, 4,3, 4,5, 5,9,
        5,10, 6,8, 7,8, 7,9, 8,11, 9,11, 9,12, 10,9}, COL);

    g.initEdges(edge_buffer);

    auto order = topSort(g);
    bool result = validateTopologicalOrder(g, order);
    if(result == false)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testGraph()
{
    testDijkstra1();
    testTopologicalSort();
}
#else
inline void testGraph(){}
#endif



#endif // _TEST_GRAPH_H_
