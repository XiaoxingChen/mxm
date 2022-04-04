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
    g.initProperty(weight_buffer);
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
    g.initProperty(weight_buffer);
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


struct MDPEdgePropertyType{
    float probability = 1.f;
    float reward = 0.f;
    bool invalid = false;
    static MDPEdgePropertyType invalidInstance()
    {
        MDPEdgePropertyType ret;
        ret.invalid = true;
        return ret;
    }
};
class MDPGraph:
    public GraphBase,
    public BinaryEdge<true, MDPGraph>,
    public BinaryEdgeProperty<MDPEdgePropertyType, MDPGraph>
{
public:
    using ThisType = MDPGraph;
    using EdgeDirectionType = BinaryEdge<true, ThisType>;
    using PropertyEdgeType = BinaryEdgeProperty<MDPEdgePropertyType, ThisType>;

    // MDPGraph(): GraphBase()
    // {
    //     this->setInvalidProperty( MDPEdgePropertyType::invalidInstance() );
    // }

    MDPGraph(size_t vertex_num): GraphBase(vertex_num)
    {
        this->setInvalidProperty( MDPEdgePropertyType::invalidInstance() );
        best_action_.resize(this->vertexNum());
    }

    MDPEdgePropertyType edgeProperty(size_t from, size_t to) const
    {
        return PropertyEdgeType::property(from, to);
    }

    ThisType& setStateVertices(const std::vector<size_t>& v)
    {
        state_vertices_ = v;
        best_action_.resize(this->vertexNum());
        return *this;
    }

    ThisType& setDiscount(float gamma)
    {
        discount_ = gamma;
        return *this;
    }

    void solve(size_t max_it)
    {
        vertex_value_ = Vector<float>::zeros(this->vertexNum());
        Vector<float> prev_vertex_value = vertex_value_;
        prev_vertex_value += 100;
        float error = 0;
        // size_t max_it = 10;
        // std::cout << "solve!!! "<< std::endl;
        for(size_t i = 0; i < max_it && !isZero(prev_vertex_value - vertex_value_, &error, 0.1); i++)
        {
            prev_vertex_value = vertex_value_;
            update();
            // std::cout << "error: " << error << std::endl;
            // std::cout << "actions: " << mxm::to_string(best_action_) << std::endl;
            // std::cout << "value: " << mxm::to_string(vertex_value_.T()) << std::endl;

        }


    }

    void update()
    {
        Vector<float> prev_vertex_value = vertex_value_;
        for(auto & state_idx: state_vertices_)
        {
            if(this->adjacency_lists.at(state_idx).empty()) continue;
            size_t best_action = this->adjacency_lists.at(state_idx).front();
            float max_action_reward = 0;
            for(auto & action_idx: this->adjacency_lists.at(state_idx))
            {
                float expected_action_reward = 0.f;
                for(auto & next_state_idx: this->adjacency_lists.at(action_idx))
                {
                    float local_reward = prev_vertex_value(next_state_idx) * discount_ + edgeProperty(action_idx, next_state_idx).reward;
                    expected_action_reward += edgeProperty(action_idx, next_state_idx).probability * local_reward;
                }
                if(expected_action_reward > max_action_reward)
                {
                    max_action_reward = expected_action_reward;
                    best_action = action_idx;
                }
            }
            vertex_value_(state_idx) = max_action_reward;
            best_action_.at(state_idx) = best_action;
        }
    }

    const Vector<float>& value() const { return vertex_value_; }

protected:
    std::vector<size_t> state_vertices_;
    std::vector<size_t> best_action_;
    Vector<float> vertex_value_;
    float discount_ = 1.f;
};

// example:
// https://youtu.be/i0o-ui1N35U?t=3717
inline void testMarkovDecisionProcess01()
{

    MDPGraph g(7);
    // 0: cool
    // 1: warm
    // 2: overheated
    // 3: cool -> slow
    // 4: cool -> fast
    // 5: warm -> slow
    // 6: warm -> fast
    Matrix<size_t> edges(fixRow(2), {
        0,3, 0,4, 1,5, 1,6,
        3,0, 4,0, 4,1,
        5,0, 5,1, 6,2}, COL);

    Vector<MDPEdgePropertyType> edge_info
    {
        {1.f, 0.f}, {1.f, 0.f}, {1.f, 0.f}, {1.f, 0.f},
        {1.f, 1.f}, {.5f, 2.f}, {.5f, 2.f},
        {.5f, 1.f}, {.5f, 1.f}, {1.f -10.f}
    };

    g.initEdges(edges);
    g.initProperty(edge_info);
    g.setStateVertices({0,1,2}).setDiscount(1);
    g.solve(2);
    if(abs(g.value()(0) - 3.5f) > eps<float>())
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
    if(abs(g.value()(1) - 2.5f) > eps<float>())
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}


inline void testFlowNetwork01()
{
    FlowNetwork<float> g(3);
    Matrix<size_t> edges(fixRow(2), {0,1, 1,2, 2,3}, COL);
    Vector<CapacityFlow<float>> properties{{3,1}, {5,1}, {4,1}};
    g.initEdges(edges);
    g.initProperty(properties);
    g.setSource(0).setSink(3);
    if(g.property(0,1).flow != 1)
    {
        std::cout << g.property(0,1).flow << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

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
    testFlowNetwork01();
    testMarkovDecisionProcess01();
}
#else
inline void testGraph(){}
#endif



#endif // _TEST_GRAPH_H_
