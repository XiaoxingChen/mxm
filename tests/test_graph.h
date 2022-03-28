#if !defined(_TEST_GRAPH_H_)
#define _TEST_GRAPH_H_

#include "test_config.h"

#if TEST_AVAILABLE_ALL
#include <iostream>
#include "mxm/graph_dijkstra.h"
#include "mxm/graph_top_sort.h"
#include "mxm/graph_utils.h"
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
    float distance = 0;

    // auto best_pred =
    std::vector<size_t> path = dijkstra(g, 1, 5, &distance);
    // auto distance = pathFromBestPredecessor(g, best_pred, 5, &path);

    if(path != std::vector<size_t>{1,3,6,5})
    {
        // std::cout << "best_pred: " << mxm::to_string(best_pred) << std::endl;
        std::cout << "path: " << mxm::to_string(path) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(distance != 20)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testDijkstra2()
{
    WeightedDirectedGraph<float> g(4);
    Matrix<size_t> edge_buffer(fixRow(2), {0,1, 1,2, 2,3, 3,0, 0,3, 0,0}, COL);
    Vector<float> weight_buffer{1,1,1,1,1,1};
    g.initEdges(edge_buffer);
    g.initWeight(weight_buffer);
    float distance = 0;

    std::vector<size_t> path = dijkstra(g, 0, 3, &distance);
    if(path != std::vector<size_t>{0,3})
    {
        // std::cout << "best_pred: " << mxm::to_string(best_pred) << std::endl;
        std::cout << "path: " << mxm::to_string(path) << std::endl;
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
        for(const auto& succ: g.adjacency(target))
        {
            if(built.count(succ) > 0)
            {
                std::cout << "target: " << target << ", succ: " << succ << std::endl;
                return false;
            }
        }
        built.insert(target);
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

    auto kahn_order = kahnAlgorithm(g);
    result == validateTopologicalOrder(g, kahn_order);
    if(result == false)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testTopologicalSort02()
{
    WeightedDirectedGraph<float> g(4);
    Matrix<size_t> edge_buffer(fixRow(2), {
        0,1, 1,2, 2,3, 3,0}, COL);

    g.initEdges(edge_buffer);

    auto order = topSort(g);
    // std::cout << mxm::to_string(order) << std::endl;
    bool result = validateTopologicalOrder(g, order);
    if(result == true)
    {
        std::cout << "this graph is expected to be unordered" << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    auto kahn_order = kahnAlgorithm(g);
    if(kahn_order.size() > 0)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

// test case source:
// https://youtu.be/7fujbpJ0LB4?list=PLDV1Zeh2NRsDGO4--qE8yH72HFL1Km93P&t=370
inline void testConnectedComponents01()
{
    UndirectedGraph g(18);
    Matrix<size_t> edges(fixRow(2), {
        6,7, 6,11, 7,11,
        4,8, 8,0, 0,4, 8,14, 0,14, 13,14, 0,13,
        5,1, 5,16, 5,17,
        15,9, 15,2, 15,10, 3,9, 2,9}, COL);

    g.initEdges(edges);
    auto result = connectedComponents(g);
    if(result.size() != 5)
    {
        std::cout << "result size: " << result.size() << std::endl;
        std::cout << mxm::to_string(result) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testWeaklyConnectedComponents01()
{
    UndirectedGraph g(7);
    Matrix<size_t> edges(fixRow(2), {
        4,0, 4,1, 5,2, 5,3, 6,4, 6,5}, COL);

    g.initEdges(edges);
    auto result = connectedComponents(g);
    if(result.size() != 1)
    {
        std::cout << "result size: " << result.size() << std::endl;
        std::cout << mxm::to_string(result) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testWeaklyConnectedComponents02()
{
    UndirectedGraph g(14);
    Matrix<size_t> edges(fixRow(2), {
        4,0, 4,1, 5,2, 5,3, 6,4, 6,5,
        11,7, 11,8, 12,9, 12,10, 13,11, 13,12}, COL);

    g.initEdges(edges);
    auto result = connectedComponents(g);
    if(result.size() != 2)
    {
        std::cout << "result size: " << result.size() << std::endl;
        std::cout << mxm::to_string(result) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testCyclic()
{
    UndirectedGraph g(14);
    Matrix<size_t> edges(fixRow(2), {
        4,0, 4,1, 5,2, 5,3, 6,4, 6,5,
        11,7, 11,8, 12,9, 12,10, 13,11, 13,12}, COL);

    g.initEdges(edges);
    auto result = isCyclic(g, connectedComponents(g).size());
    if(result == true)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testCyclic01()
{

}

inline void testGraph()
{
    testDijkstra1();
    testDijkstra2();
    testTopologicalSort();
    testTopologicalSort02();
    testConnectedComponents01();
    testWeaklyConnectedComponents01();
    testWeaklyConnectedComponents02();
}
#else
inline void testGraph(){}
#endif



#endif // _TEST_GRAPH_H_
