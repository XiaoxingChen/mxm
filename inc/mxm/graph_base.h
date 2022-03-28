#if !defined(_GRAPH_BASE_H_)
#define _GRAPH_BASE_H_

#include "linalg_mat.h"
#include "common.h"
#include <map>
#include <type_traits>

namespace mxm
{

class GraphBase
{
public:
    using ThisType = GraphBase;
    const size_t& vertexNum() const { return vertex_num_; }
    GraphBase() {}
    GraphBase(size_t vertex_num): vertex_num_(vertex_num)  {}
    ThisType& setEdgeBuffer(const Matrix<size_t>& edges)
    {
        edge_buffer_ = edges;
    }
    ThisType& addEdge(const Vector<size_t>& edges)
    {
        edge_buffer_ = hstack(edge_buffer_, edges);
    }
    bool validVertex(size_t idx) const { return idx < vertex_num_; }
protected:
    size_t vertex_num_ = 0;
    Matrix<size_t> edge_buffer_;
};

template<bool Directed, typename GraphType>
class BinaryEdge
{
public:
    using ThisType = BinaryEdge<Directed, GraphType>;

    const std::vector<size_t>& adjacency(size_t v_index) const { return adjacency_lists.at(v_index); }
    // std::vector<size_t>& adjacency(size_t v_index) { return adjacency_lists.at(v_index); }

    static constexpr bool directed() { return Directed; }

    void initEdges(const Matrix<size_t>& edges)
    {
        size_t vertex_num = static_cast<GraphType*>(this)->vertexNum();
        initEdgeIndices(edges);
        initAdjacencyLists(vertex_num, edges);
    }

    // binary edge interface
    template<bool D=Directed, std::enable_if_t<D, int> T=0>
    size_t edgeIndex(size_t from, size_t to) const
    {
        if(0 == edge_index_.count({from, to})) return edge_num_ + 1;
        return edge_index_.at({from, to});
    }

    template<bool D=Directed, std::enable_if_t<!D, int> T=0>
    size_t edgeIndex(size_t v1, size_t v2) const
    {
        size_t from = std::min(v1, v2);
        size_t to = std::max(v1, v2);
        if(0 == edge_index_.count({from, to})) return edge_num_ + 1;
        return edge_index_.at({from, to});
    }

    // binary edge interface
    Vector<size_t> edgeIndices(size_t from_idx, const std::vector<size_t>& vertices) const
    {
        Vector<size_t> ret(vertices.size());
        for(size_t query_i = 0; query_i < vertices.size(); query_i++)
        {
            ret(query_i) = edgeIndex(from_idx, vertices.at(query_i));
        }
        return ret;
    }

    // multi-edge interface
    Vector<size_t> edgeIndices(const Matrix<size_t>& query_edges) const
    {
        Vector<size_t> ret(query_edges.shape(1));
        for(size_t query_i = 0; query_i < query_edges.shape(1); query_i++)
        {
            size_t from_idx = query_edges(0, query_i);
            size_t to_idx = query_edges(1, query_i);
            ret(query_i) = edgeIndex(from_idx, to_idx);
        }
        return ret;
    }

    size_t edgeNum() const { return edge_num_; }

protected:
    size_t edge_num_ = 0;
    std::vector<std::vector<size_t>> adjacency_lists;
    std::map<std::array<size_t, 2>, size_t > edge_index_;

private:

    template<bool D=Directed, std::enable_if_t<D, int> T=0>
    void initAdjacencyLists(size_t vertex_num, const Matrix<size_t>& edge_buffer)
    {
        adjacency_lists.resize(vertex_num);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            adjacency_lists.at(edge_buffer(0, i)).push_back(edge_buffer(1, i));
        }
    }

    template<bool D=Directed, std::enable_if_t<!D, int> T=0>
    void initAdjacencyLists(size_t vertex_num, const Matrix<size_t>& edge_buffer)
    {
        adjacency_lists.resize(vertex_num);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            adjacency_lists.at(edge_buffer(0, i)).push_back(edge_buffer(1, i));
            adjacency_lists.at(edge_buffer(1, i)).push_back(edge_buffer(0, i));
        }
    }

    template<bool D=Directed, std::enable_if_t<D, int> T=0>
    void initEdgeIndices(const Matrix<size_t>& edge_buffer)
    {
        edge_num_ = edge_buffer.shape(1);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            edge_index_[{edge_buffer(0, i), edge_buffer(1, i)}] = i;
        }
    }

    template<bool D=Directed, std::enable_if_t<!D, int> T=0>
    void initEdgeIndices(const Matrix<size_t>& edge_buffer)
    {
        edge_num_ = edge_buffer.shape(1);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            size_t from = std::min(edge_buffer(0, i), edge_buffer(1, i));
            size_t to = std::max(edge_buffer(0, i), edge_buffer(1, i));
            edge_index_[{from, to}] = i;
        }
    }


};

template<typename WeightType>
class WeightedEdge
{
public:
    void setInvalidWeight(WeightType invalid_weight)
    {
        invalid_weight_ = invalid_weight;
    }

    void initWeight(const Vector<WeightType>& w_buffer)
    {
        weight_buffer_ = w_buffer;
    }

    WeightType weight(size_t edge_index) const
    {
        if(edge_index >= weight_buffer_.size()) return invalid_weight_;
        return weight_buffer_(edge_index);
    }
protected:
    Vector<WeightType> weight_buffer_;
    WeightType invalid_weight_;
};

#if 0
#else


class DirectedGraph: public GraphBase, public BinaryEdge<true, DirectedGraph>
{
public:
    DirectedGraph(): GraphBase()
    {
    }

    DirectedGraph(size_t vertex_num): GraphBase(vertex_num)
    {
    }

private:
};

template<typename DType>
class WeightedDirectedGraph:
    public GraphBase,
    public BinaryEdge<true, WeightedDirectedGraph<DType>>,
    public WeightedEdge<DType>
{
public:
    using DistanceType = DType;
    using EdgeDirectionType = BinaryEdge<true, WeightedDirectedGraph<DType>>;

    WeightedDirectedGraph(): GraphBase()
    {
        this->setInvalidWeight(INFINITY);
    }

    WeightedDirectedGraph(size_t vertex_num): GraphBase(vertex_num)
    {
        this->setInvalidWeight(INFINITY);
    }

    DType weight(size_t from, size_t to) const
    {
        size_t edge_index = this->edgeIndex(from, to);
        return WeightedEdge<DType>::weight(edge_index);
    }
protected:
};

class UndirectedGraph:
    public GraphBase,
    public BinaryEdge<false, UndirectedGraph>
{
public:
    UndirectedGraph(): GraphBase()
    {
    }

    UndirectedGraph(size_t vertex_num): GraphBase(vertex_num)
    {
    }

protected:
};

#endif


// template<typename ExtendedType>
// class IndexExtension
// {
// public:
//     size_t rawFromExt(const ExtendedType& ext);
//     ExtendedType extFromRaw(size_t raw);
// };

class GridGraph
{
public:
private:
};

template<class T, class = void>
struct is_weighted_binary_edge: std::false_type{};

template<class T>
struct is_weighted_binary_edge<
    T, typename mxm::void_t< decltype(std::declval<T>().weight(0,0)) >::type
    >: std::true_type{};


} // namespace mxm


#endif // _GRAPH_BASE_H_
