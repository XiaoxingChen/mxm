#if !defined(_TEST_GRAPH_H_)
#define _TEST_GRAPH_H_

#include "test_config.h"

#if TEST_AVAILABLE_ALL
#include <iostream>
#include "mxm/graph_dijkstra.h"
#include "mxm/graph_top_sort.h"
#include "mxm/graph_utils.h"
#include "mxm/random.h"
#include "mxm/graph_flow.h"
#endif

using namespace mxm;
#if TEST_AVAILABLE_ALL
// test graph:
// link: https://en.wikipedia.org/wiki/File:Dijkstra_Animation.gif
inline void testShortedPath1()
{
    WeightedDirectedGraph<float> g(7);
    Matrix<size_t> edge_buffer(fixRow(2), {1,2, 1,3, 1,6, 2,4, 2,3, 3,6, 3,4, 6,5, 4,5}, COL);
    Vector<float> weight_buffer{7,9,14,15,10,2,11,9,6};
    g.initEdges(edge_buffer);
    g.initProperty(weight_buffer);
    float distance = 0;

    for(size_t method_idx = 0; method_idx < 2; method_idx++)
    {
        std::vector<size_t> path;
        if(0 == method_idx) path = dijkstra(g, 1, 5, &distance);
        else if(1 == method_idx) path = shortedPathBellmanFord(g, 1, 5, &distance);

        if(path != std::vector<size_t>{1,3,6,5})
        {
            // std::cout << "best_pred: " << mxm::to_string(best_pred) << std::endl;
            std::cout << "method_idx: " << method_idx << ", path: " << mxm::to_string(path) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(distance != 20)
        {
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

}

inline void testShortedPath2()
{
    WeightedDirectedGraph<float> g(4);
    Matrix<size_t> edge_buffer(fixRow(2), {0,1, 1,2, 2,3, 3,0, 0,3, 0,0}, COL);
    Vector<float> weight_buffer{1,1,1,1,1,1};
    g.initEdges(edge_buffer);
    g.initProperty(weight_buffer);
    float distance = 0;

    for(size_t method_idx = 0; method_idx < 2; method_idx++)
    {
        std::vector<size_t> path;
        if(0 == method_idx) path = dijkstra(g, 0, 3, &distance);
        else if(1 == method_idx) path = shortedPathBellmanFord(g, 0, 3, &distance);

        if(path != std::vector<size_t>{0,3})
        {
            // std::cout << "best_pred: " << mxm::to_string(best_pred) << std::endl;
            std::cout << "method_idx: " << method_idx << ", path: " << mxm::to_string(path) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }


}

// test graph:
// link: https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/
inline void testShortedPath3()
{
    WeightedUnirectedGraph<float> g(9);
    Matrix<size_t> edge_buffer(fixRow(2), {0,1, 0,7, 1,2, 1,7, 2,3, 2,5, 2,8, 3,5, 3,4, 4,5, 5,6, 6,7, 6,8, 7,8}, COL);
    Vector<float> weight_buffer{4, 8, 8, 11, 7, 4, 2, 9, 14, 10, 2, 1, 6, 7};
    g.initEdges(edge_buffer);
    g.initProperty(weight_buffer);
    float distance = 0;

    for(size_t method_idx = 0; method_idx < 2; method_idx++)
    {
        std::vector<size_t> path;
        if(0 == method_idx) path = dijkstra(g, 0, 4, &distance);
        else if(1 == method_idx) path = shortedPathBellmanFord(g, 0, 4, &distance);

        if(path != std::vector<size_t>{0,7,6,5,4})
        {
            std::cout << "distance: " << distance << std::endl;
            std::cout << "path: " << mxm::to_string(path) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

inline void testShortedPath4()
{
    WeightedUnirectedGraph<float> g(4);
    Matrix<size_t> edge_buffer(fixRow(2), {0,1, 1,2, 1,3, 2,3}, COL);
    Vector<float> weight_buffer{1, 2, 10, 9};
    g.initEdges(edge_buffer);
    g.initProperty(weight_buffer);
    float distance = 0;

    for(size_t method_idx = 0; method_idx < 2; method_idx++)
    {
        std::vector<size_t> path;
        if(0 == method_idx) path = dijkstra(g, 0, 3, &distance);
        else if(1 == method_idx) path = shortedPathBellmanFord(g, 0, 3, &distance);

        // std::vector<size_t> path = dijkstra(g, 0, 3, &distance);
        if(path != std::vector<size_t>{0,1,3})
        {
            std::cout << "distance: " << distance << std::endl;
            std::cout << "path: " << mxm::to_string(path) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

inline void testShortedPath5()
{
    WeightedUnirectedGraph<float> g(4);
    Matrix<size_t> edge_buffer(fixRow(2), {0,1, 2,3}, COL);
    Vector<float> weight_buffer{1, 2};
    g.initEdges(edge_buffer);
    g.initProperty(weight_buffer);
    float distance = 0;

    for(size_t method_idx = 0; method_idx < 2; method_idx++)
    {
        std::vector<size_t> path;
        if(0 == method_idx) path = dijkstra(g, 0, 3, &distance);
        else if(1 == method_idx) path = shortedPathBellmanFord(g, 0, 3, &distance);

        // std::vector<size_t> path = dijkstra(g, 0, 3, &distance);
        if(path != std::vector<size_t>{})
        {
            std::cout << "distance: " << distance << std::endl;
            std::cout << "path: " << mxm::to_string(path) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
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

struct MDPTransition{
    size_t next_state_idx = 0;
    float probability = 1.f;
    float reward = 0.f;
};
class MarkovDecisionProcessSolver
{
public:
    using ThisType = MarkovDecisionProcessSolver;

    ThisType& setQValueTable(const Matrix<std::vector<MDPTransition>>& table)
    {
        qValueTable_ = table;
        state_value_ = Vector<float>::zeros(stateNum());
        state_best_action_ = Vector<size_t>::zeros(stateNum());
        return *this;
    }
    size_t stateNum() const { return qValueTable_.shape(1); }
    size_t actionNum() const { return qValueTable_.shape(0); }

    ThisType& setDiscount(float gamma)
    {
        discount_ = gamma;
        return *this;
    }

    void qLearningSolve(size_t init_state_idx, size_t exit_state_idx, size_t max_it=100)
    {
        size_t state_idx = init_state_idx;
        size_t prev_state_idx = init_state_idx;
        size_t prev_action_idx = 0;
        size_t prev_transition_idx = 0;

        std::vector<std::string> node_strs = {"(1,A)", "(2,A)", "(3,A)", "(1,B)", "(2,B)", "(3,B)", "(1,C)", "(2,C)", "(3,C)", "(1,D)", "(2,D)", "(3,D)"};
        std::vector<std::string> action_strs = {"v", "^", ">", "<"};

        for(size_t it_cnt = 0; it_cnt < max_it; it_cnt++)
        {
            if(state_idx == exit_state_idx) break;
            float max_action_value = 0.f;
            size_t best_action = 0;
            Vector<float> prov_vec;

            for(size_t action_idx = 0; action_idx < actionNum(); action_idx++)
            {
                float action_expect_value = 0;
                std::vector<float> prob_vec_data;
                for(auto & transition : qValueTable_(action_idx, state_idx))
                {
                    action_expect_value += transition.probability*(state_value_(transition.next_state_idx) * discount_ + transition.reward);
                    prob_vec_data.push_back(transition.probability);
                }
                if(action_expect_value > max_action_value)
                {
                    best_action = action_idx;
                    max_action_value = action_expect_value;
                    prov_vec = Vector<float>(std::move(prob_vec_data));
                }
            }
            if(max_action_value < eps<float>())
            {
                best_action = random::weightedSample<float>(Vector<float>::ones(actionNum()) / actionNum());
            }

            size_t transition_idx = random::weightedSample<float>(prov_vec);
            float new_q_s_a_val = 0;

            if( it_cnt > 0 )
            {
                qValueTable_(prev_action_idx, prev_state_idx).at(prev_transition_idx).reward = state_value_(state_idx) + discount_ * qValueTable_(best_action, state_idx).at(transition_idx).reward;
                new_q_s_a_val = qValueTable_(prev_action_idx, prev_state_idx).at(prev_transition_idx).reward;
            }


            prev_action_idx = best_action;
            prev_state_idx = state_idx;
            prev_transition_idx = transition_idx;
            state_idx = qValueTable_(best_action, state_idx).at(transition_idx).next_state_idx;
            std::cout << "i: " << it_cnt << ", action: " << action_strs.at(best_action) << ", subsequent state: " << node_strs.at(state_idx) << ", Q" << node_strs.at(prev_state_idx) << action_strs.at(prev_action_idx) << ": " << new_q_s_a_val << std::endl;


        }
    }

    void solve(size_t max_it)
    {
        state_value_ *= 0.f;
        Vector<float> prev_state_value = state_value_;
        prev_state_value += 100;
        float error = 0;
        // size_t max_it = 10;
        // std::cout << "solve!!! "<< std::endl;
        for(size_t i = 0; i < max_it && !isZero(prev_state_value - state_value_, &error, 0.1); i++)
        {
            prev_state_value = state_value_;
            update();
            // std::cout << "error: " << error << std::endl;
            // std::cout << "actions: " << mxm::to_string(state_best_action_.T()) << std::endl;
            // std::cout << "value: " << mxm::to_string(state_value_.T()) << std::endl;

        }

    }
    void update()
    {
        Vector<float> shadow_state_value = state_value_;
        for(size_t state_idx = 0; state_idx < stateNum(); state_idx++)
        {
            float max_action_value = 0.f;
            size_t best_action = 0;
            for(size_t action_idx = 0; action_idx < actionNum(); action_idx++)
            {
                float action_expect_value = 0;
                for(auto & transition : qValueTable_(action_idx, state_idx))
                {
                    action_expect_value += transition.probability*(shadow_state_value(transition.next_state_idx) * discount_ + transition.reward);
                }
                if(action_expect_value > max_action_value)
                {
                    best_action = action_idx;
                    max_action_value = action_expect_value;
                }
            }
            state_value_(state_idx) = max_action_value;
            state_best_action_(state_idx) = best_action;
        }
    }

    ThisType& setStateValue(const Vector<float> state_value)
    {
        state_value_ = state_value_;
        return *this;
    }
private:
    Matrix<std::vector<MDPTransition>> qValueTable_;
    Vector<float> state_value_;
    Vector<size_t> state_best_action_;
    float discount_ = 1.f;
};


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
            if(this->adjacency_lists_.at(state_idx).empty()) continue;
            size_t best_action = this->adjacency_lists_.at(state_idx).front();
            float max_action_reward = 0;
            for(auto & action_idx: this->adjacency_lists_.at(state_idx))
            {
                float expected_action_reward = 0.f;
                for(auto & next_state_idx: this->adjacency_lists_.at(action_idx))
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

// example:
// https://youtu.be/i0o-ui1N35U?t=3717
inline void testMarkovDecisionProcess02()
{
    MarkovDecisionProcessSolver mdp_solver;
    Matrix<std::vector<MDPTransition>> qValueTable({2,3},
    {
        {/*cool slow*/{0, 1.f, 1.f}}, {/*warm slow*/{0, 0.5, 1}, {1, 0.5, 1}},{}, //slow
        {/*cool fast*/{0, 0.5f, 2.f},{1, 0.5f, 2.f}},{/*warm fast*/{2, 1.f, -10.f}},{} //fast
    });
    mdp_solver.setDiscount(1.f);
    mdp_solver.setQValueTable(qValueTable);
    mdp_solver.solve(3);
}

inline void testMarkovDecisionProcess03()
{
    MarkovDecisionProcessSolver mdp_solver;
    Matrix<std::vector<MDPTransition>> qValueTable({4,12},
    {
        // 1A               2A               [2]:3A           [3]:1B          [4]:2B         [5]:3B           [6]:1C        [7]:2C          [8]:3C           [9]:1D          [10]:2D          [11]:3D
        {{1, 1.f, .2f}},{{2, 1.f, .2f}},{{2, 1.f, .0f}},{{3, 1.f, .0f}},{{5, 1.f, .0f}},{{5, 1.f, 0.f}},{{7, 1.f, .1f}},{{8, 1.f, 0.f}},{{8, 1.f, 0.f}},{{10,1.f, .1f}},{{11, 1.f, .1f}},{{11, 1.f, 0.f}},
        {{0, 1.f, .0f}},{{1, 1.f, .0f}},{{2, 1.f, .0f}},{{4, 1.f, .0f}},{{3, 1.f, .1f}},{{4, 1.f, .1f}},{{6, 1.f, 0.f}},{{6, 1.f, 0.f}},{{7, 1.f, 0.f}},{{9, 1.f, 0.f}},{{9,  1.f, 0.f}},{{10, 1.f, 0.f}},
        {{0, 1.f, .0f}},{{1, 1.f, .0f}},{{5, 1.f, .2f}},{{6, 1.f, .1f}},{{7, 1.f, .2f}},{{8, 1.f, .2f}},{{9, 1.f, 0.f}},{{10,1.f, .2f}},{{5, 1.f, 0.f}},{{9, 1.f, 0.f}},{{10, 1.f, 0.f}},{{11, 1.f, 0.f}},
        {{0, 1.f, .0f}},{{1, 1.f, .0f}},{{2, 1.f, .0f}},{{3, 1.f, .0f}},{{4, 1.f, .0f}},{{2, 1.f, 0.f}},{{3, 1.f, .1f}},{{4, 1.f, 0.f}},{{8, 1.f, 0.f}},{{6, 1.f, 0.f}},{{7,  1.f, 0.f}},{{11, 1.f, 0.f}}
    }, ROW);
    Vector<float> state_value = Vector<float>::zeros(12);
    state_value(11) = 1;
    mdp_solver.setStateValue(state_value);
    mdp_solver.setDiscount(.9f);
    mdp_solver.setQValueTable(qValueTable);
    // mdp_solver.qLearningSolve(0, 11, 20);
}


inline void testFlowNetwork01()
{
    FlowNetwork<float> g(4);
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

inline void testFlowNetwork02()
{
    WeightedDirectedGraph<float> g(4);
    Matrix<size_t> edge_buffer(fixRow(2), {
        0,1, 0,2, 2,1, 1,3, 2,3}, COL);

    g.initEdges(edge_buffer);
    Vector<float> cap{10, 3, 3, 10, 2};
    g.initProperty(cap);

    fordFulkersonMaxFLow(g, 0, 3);


}

inline void testGraph()
{
    testShortedPath1();
    testShortedPath2();
    testShortedPath3();
    testShortedPath4();
    testShortedPath5();
    testTopologicalSort();
    testTopologicalSort02();
    testConnectedComponents01();
    testWeaklyConnectedComponents01();
    testWeaklyConnectedComponents02();
    testFlowNetwork01();
    testMarkovDecisionProcess01();
    testMarkovDecisionProcess02();
    testMarkovDecisionProcess03();
    testFlowNetwork02();
}
#else
inline void testGraph(){}
#endif



#endif // _TEST_GRAPH_H_
