#if !defined(_TEST_GRAPH_H_)
#define _TEST_GRAPH_H_

#include "test_config.h"

#if TEST_AVAILABLE_ALL
#include <iostream>
#include "mxm/graph_dijkstra.h"
#endif

using namespace mxm;
#if TEST_AVAILABLE_ALL
// test graph:
// link: https://en.wikipedia.org/wiki/File:Dijkstra_Animation.gif
inline void testDijkstra1()
{
    SparseDirectedWeightedGraph g(7);
    g.addEdge(1,2,7);
    g.addEdge(1,3,9);
    g.addEdge(1,6,14);
    g.addEdge(2,4,15);
    g.addEdge(2,3,10);
    g.addEdge(3,6,2);
    g.addEdge(3,4,11);
    g.addEdge(6,5,9);
    g.addEdge(4,5,6);

    auto best_pred = dijkstra(g, 1);
    std::vector<size_t> path;
    auto distance = pathFromBestPredecessor(g, best_pred, 5, path);
    if(distance != 20)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
    if(path != std::vector<size_t>{1,3,6,5})
    {
        std::cout << mxm::to_string(path) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testGraph()
{
    testDijkstra1();
}
#else
inline void testGraph(){}
#endif



#endif // _TEST_GRAPH_H_
