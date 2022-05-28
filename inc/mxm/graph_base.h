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

    bool validVertex(size_t idx) const { return idx < vertex_num_; }
    static constexpr size_t nullVertex() { return std::numeric_limits<size_t>::max(); }
protected:
    size_t vertex_num_ = 0;
};

template<bool Directed, typename GraphType>
class BinaryEdge
{
public:
    using ThisType = BinaryEdge<Directed, GraphType>;

    const std::vector<size_t>& adjacency(size_t v_index) const { return adjacency_lists_.at(v_index); }
    // std::vector<size_t>& adjacency(size_t v_index) { return adjacency_lists_.at(v_index); }

    static constexpr bool directed() { return Directed; }

    void initEdges(const Matrix<size_t>& edges)
    {
        size_t vertex_num = static_cast<GraphType*>(this)->vertexNum();
        initEdgeIndices(edges);
        initAdjacencyLists(vertex_num, edges);
    }

    const std::map<std::array<size_t, 2>, size_t >&
    edgeIndices() const { return edge_index_; }

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
    std::vector<std::vector<size_t>> adjacency_lists_;
    std::map<std::array<size_t, 2>, size_t > edge_index_;

private:

    template<bool D=Directed, std::enable_if_t<D, int> T=0>
    void initAdjacencyLists(size_t vertex_num, const Matrix<size_t>& edge_buffer)
    {
        adjacency_lists_.resize(vertex_num);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            adjacency_lists_.at(edge_buffer(0, i)).push_back(edge_buffer(1, i));
        }
    }

    template<bool D=Directed, std::enable_if_t<!D, int> T=0>
    void initAdjacencyLists(size_t vertex_num, const Matrix<size_t>& edge_buffer)
    {
        adjacency_lists_.resize(vertex_num);
        for(size_t i = 0; i < edge_buffer.shape(1); i++)
        {
            adjacency_lists_.at(edge_buffer(0, i)).push_back(edge_buffer(1, i));
            adjacency_lists_.at(edge_buffer(1, i)).push_back(edge_buffer(0, i));
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

template<bool IsDirected>
class Graph:
    public GraphBase,
    public BinaryEdge<IsDirected, Graph<IsDirected>>
{
public:
    Graph(): GraphBase()
    {
    }

    Graph(size_t vertex_num): GraphBase(vertex_num)
    {
    }

protected:
};

using UndirectedGraph = Graph<false>;
using DirectedGraph = Graph<true>;

template<typename PropertyType, bool IsDirected>
class EdgeProperty
{
public:
    using ThisType = EdgeProperty<PropertyType, IsDirected>;
    ThisType& setInvalidProperty(PropertyType invalid_property)
    {
        invalid_property_ = invalid_property;
        return *this;
    }

    ThisType& initProperty(const Vector<PropertyType>& property_buffer)
    {
        property_buffer_ = property_buffer;
        return *this;
    }

    ThisType& setTopology( const Graph<IsDirected>* p )
    {
        p_topology_ = p;
        return *this;
    }

    const Graph<IsDirected>*
    topology() const { return p_topology_; }

    const PropertyType& property(size_t edge_index) const
    {
        assert(edge_index < property_buffer_.size());
        if(edge_index >= property_buffer_.size()) return invalid_property_;
        return property_buffer_(edge_index);
    }

    PropertyType& property(size_t edge_index)
    {
        assert(edge_index < property_buffer_.size());
        if(edge_index >= property_buffer_.size())
            return invalid_property_;

        return property_buffer_(edge_index);
    }

    const Vector<PropertyType>& properties() const
    {
        return property_buffer_;
    }

protected:
    Vector<PropertyType> property_buffer_;
    PropertyType invalid_property_;
    const Graph<IsDirected>* p_topology_ = nullptr;
};

template<typename PropertyType, bool IsDirected>
class BinaryEdgeProperty: public EdgeProperty<PropertyType, IsDirected>
{
public:
    const PropertyType& operator()(size_t v1, size_t v2) const
    {
        assert(this->p_topology_);
        size_t edge_index = this->p_topology_->edgeIndices(Vector<size_t>{v1, v2})(0);
        return EdgeProperty<PropertyType, IsDirected>::property(edge_index);
    }

    PropertyType& operator()(size_t v1, size_t v2)
    {
        assert(this->p_topology_);
        size_t edge_index = this->p_topology_->edgeIndices(Vector<size_t>{v1, v2})(0);
        return EdgeProperty<PropertyType, IsDirected>::property(edge_index);
    }

protected:

};


} // namespace mxm


#endif // _GRAPH_BASE_H_
