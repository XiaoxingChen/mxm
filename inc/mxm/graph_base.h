#if !defined(_GRAPH_BASE_H_)
#define _GRAPH_BASE_H_

#include "linalg_mat.h"
#include <map>

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

class DirectedBinaryEdge
{
public:
    using ThisType = DirectedBinaryEdge;

    const std::vector<size_t>& successors(size_t v_index) const { return successors_.at(v_index); }
    std::vector<size_t>& successors(size_t v_index) { return successors_.at(v_index); }

    void initSuccessorList(size_t vertex_num, const Matrix<size_t>& edge_buffer)
    {
        successors_.resize(vertex_num);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            successors_.at(edge_buffer(0, i)).push_back(edge_buffer(1, i));
        }
    }

    void initEdgeIndices(const Matrix<size_t>& edge_buffer)
    {
        edge_num_ = edge_buffer.shape(1);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            edge_index_[{edge_buffer(0, i), edge_buffer(1, i)}] = i;
        }
    }

    // binary edge interface
    size_t edgeIndex(size_t from, size_t to) const
    {
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

protected:
    size_t edge_num_ = 0;
    std::vector<std::vector<size_t>> successors_;
    std::map<std::array<size_t, 2>, size_t > edge_index_;

};


class UndirectedBinaryEdge
{
public:
    using ThisType = UndirectedBinaryEdge;

    const std::vector<size_t>& neighbors(size_t v_index) const { return neighbors_.at(v_index); }
    std::vector<size_t>& neighbors(size_t v_index) { return neighbors_.at(v_index); }

    void initSuccessorList(size_t vertex_num, const Matrix<size_t>& edge_buffer)
    {
        neighbors_.resize(vertex_num);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            neighbors_.at(edge_buffer(0, i)).push_back(edge_buffer(1, i));
            neighbors_.at(edge_buffer(1, i)).push_back(edge_buffer(0, i));
        }
    }

protected:
    std::vector<std::vector<size_t>> neighbors_;
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
class SparseDirectedWeightedGraph
{
public:
    using DType = FloatType;
    SparseDirectedWeightedGraph(size_t node_num)
    :successors_(node_num)
    {

    }

    void addEdge(size_t from, size_t to, DType distance)
    {
        if(from >= nodeNum()) resize(from + 1);
        successors_.at(from).insert({distance, to});
    }

    const std::multimap<DType, size_t>& successors(size_t node_idx) const
    {
        return successors_.at(node_idx);
    }

    void resize(size_t s)
    {
        successors_.resize(s);
    }

    bool valid(size_t idx) const { return idx < nodeNum(); }

    DType distance(size_t from, size_t to) const
    {
        if(from == to) return 0;
        if( ! valid(from) ) return INFINITY;
        for(auto & pair: successors_.at(from))
        {
            if(pair.second == to) return pair.first;
        }
        return INFINITY;
    }

    size_t nodeNum() const { return successors_.size(); }

private:
    std::vector<std::multimap<DType, size_t>> successors_;

};
#else


class DirectedGraph: public GraphBase, public DirectedBinaryEdge
{
public:
    DirectedGraph(): GraphBase()
    {
    }

    DirectedGraph(size_t vertex_num): GraphBase(vertex_num)
    {
    }

    void initEdges(const Matrix<size_t>& edges)
    {
        initEdgeIndices(edges);
        initSuccessorList(vertexNum(), edges);
    }
private:
};

template<typename DType>
class WeightedDirectedGraph: public GraphBase, public DirectedBinaryEdge, public WeightedEdge<DType>
{
public:
    using DistanceType = DType;

    WeightedDirectedGraph(): GraphBase()
    {
        this->setInvalidWeight(INFINITY);
    }

    WeightedDirectedGraph(size_t vertex_num): GraphBase(vertex_num)
    {
        this->setInvalidWeight(INFINITY);
    }

    void initEdges(const Matrix<size_t>& edges)
    {
        initEdgeIndices(edges);
        initSuccessorList(vertexNum(), edges);
    }

    DType distance(size_t from, size_t to) const
    {
        size_t edge_index = edgeIndex(from, to);
        return this->weight(edge_index);
    }
protected:
};

#endif


#if 0
// fast indexing
template<typename WeightType>
class BinaryEdgeGraph
{
public:
    const size_t& vertexNum() const { return vertex_num_; }
    size_t& vertexNum() { return vertex_num_; }

    GeneralGraph& setEdgeBuffer(const Matrix<size_t>& edges)
    {

    }

    GeneralGraph& setWeightBuffer(const Vector<WeightType>& weights)
    {

    }

    const WeightType& weight(const Vector<size_t>& edge)
    {

    }



private:
    size_t vertex_num_;
    std::vector<std::vector<size_t>> successors_;

    Vector<WeightType> weight_buffer_;
};
#endif
} // namespace mxm


#endif // _GRAPH_BASE_H_
