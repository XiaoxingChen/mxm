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
#if 0
    ThisType& setEdgeBuffer(const Matrix<size_t>& edges)
    {
        edge_buffer_ = edges;
        return *this;
    }
    ThisType& addEdge(const Vector<size_t>& edges)
    {
        edge_buffer_ = hstack(edge_buffer_, edges);
        return *this;
    }

    const Matrix<size_t>& edgeBuffer() const { return edge_buffer_; }
#endif

    bool validVertex(size_t idx) const { return idx < vertex_num_; }
    static constexpr size_t nullVertex() { return std::numeric_limits<size_t>::max(); }
protected:
    size_t vertex_num_ = 0;
    // Matrix<size_t> edge_buffer_;
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
class UnweightedGraph:
    public GraphBase,
    public BinaryEdge<IsDirected, UnweightedGraph<IsDirected>>
{
public:
    UnweightedGraph(): GraphBase()
    {
    }

    UnweightedGraph(size_t vertex_num): GraphBase(vertex_num)
    {
    }

protected:
};

using UndirectedGraph = UnweightedGraph<false>;
using DirectedGraph = UnweightedGraph<true>;

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

    ThisType& setTopology( const UnweightedGraph<IsDirected>* p )
    {
        p_topology_ = p;
        return *this;
    }

    const UnweightedGraph<IsDirected>*
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
    const UnweightedGraph<IsDirected>* p_topology_ = nullptr;
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
        // if(this->p_topology_ == nullptr)
        // {
        //     exit(0);
        // }
        size_t edge_index = this->p_topology_->edgeIndices(Vector<size_t>{v1, v2})(0);
        return EdgeProperty<PropertyType, IsDirected>::property(edge_index);
    }

protected:

};




#if 0
template<typename DType>
struct CapacityFlow
{
    CapacityFlow() {}
    CapacityFlow(DType _cap, DType _flow=DType(0)):cap(_cap), flow(_flow) {}
    DType flow = DType(0);
    DType cap = DType(0);
    bool invalid = false;
    static CapacityFlow invalidInstance() { CapacityFlow obj(0,0); obj.invalid = true; return obj; }
};

template<typename DType, std::enable_if_t<std::is_floating_point<DType>::value, int> T=0>
class FlowNetwork:
    public GraphBase,
    public BinaryEdge<true, FlowNetwork<DType>>,
    public BinaryEdgeProperty<CapacityFlow<DType>, FlowNetwork<DType>>
{
public:
    using ThisType = FlowNetwork<DType>;
    using EdgeDirectionType = BinaryEdge<true, ThisType>;
    using PropertyEdgeType = BinaryEdgeProperty<CapacityFlow<DType>, ThisType>;

    FlowNetwork(): GraphBase()
    {
        this->setInvalidProperty(CapacityFlow<DType>::invalidInstance());
    }

FlowNetwork(size_t vertex_num): GraphBase(vertex_num)
    {
        this->setInvalidProperty(CapacityFlow<DType>::invalidInstance());
    }
    ThisType& setSource(size_t source) { source_ = source; return *this;}
    ThisType& setSink(size_t sink) { sink_ = sink; return *this;}

    size_t source(size_t source) const { return source_; }
    size_t sink(size_t sink) const { return sink_; }
protected:
    size_t source_;
    size_t sink_;
};


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
#endif
template<class GraphType>
struct is_edge_binary:
std::integral_constant<bool, true>
{};

#if 0
template<class T, class = void>
struct is_edge_weighted: std::false_type{};

template<class T>
struct is_edge_weighted<
    T, typename std::void_t< decltype(std::declval<T>().weight(0,0)) >
    >: std::true_type{};

template<class GraphType>
struct is_weighted_binary_edge:
std::integral_constant<bool, is_edge_binary<GraphType>::value && is_edge_weighted<GraphType>::value >
{};

template<class GraphType>
struct is_unweighted_binary_edge:
std::integral_constant<bool, is_edge_binary<GraphType>::value && !is_edge_weighted<GraphType>::value >
{};

template <typename GraphType>
inline constexpr bool is_edge_weighted_v = is_edge_weighted<GraphType>::value;

template <typename GraphType>
inline constexpr bool is_edge_binary_v = is_edge_binary<GraphType>::value;

template <typename GraphType>
inline constexpr bool is_unweighted_binary_edge_v = is_unweighted_binary_edge<GraphType>::value;
#endif
#if 0
template<typename DType, bool IsDirected, std::enable_if_t<std::is_floating_point<DType>::value, int> T=0>
class WeightedGraph
{
public:
    using ThisType = WeightedGraph<DType, IsDirected>;
    using WeightType = DType;

    WeightedGraph()
    {
        edge_property_.setInvalidProperty(std::numeric_limits<DType>::max());
    }

    WeightedGraph(size_t vertex_num)
    {
        p_topology_ = std::make_shared<UnweightedGraph<IsDirected>>(vertex_num);
        edge_property_.setInvalidProperty(std::numeric_limits<DType>::max());
    }

    WeightedGraph(std::shared_ptr<UnweightedGraph<IsDirected>> p_topology)
    {
        p_topology_ = p_topology;
        edge_property_.setInvalidProperty(std::numeric_limits<DType>::max());
    }

    void setInvalidProperty(const DType& v) { edge_property_.setInvalidProperty(v); }
    void initWeight(const Vector<WeightType>& property_buffer) { edge_property_.initProperty(property_buffer); }
    void initEdges(const Matrix<size_t>& edges) const {p_topology_->initEdges(edges); }
    static constexpr bool directed() { return IsDirected; }

    const Vector<DType>& weights() const
    {
        return edge_property_.properties();
    }

    const DType& weight(size_t from, size_t to) const
    {
        return edge_property_.property(p_topology_->edgeIndex(from, to));
    }

    DType& weight(size_t from, size_t to)
    {
        return edge_property_.property(p_topology_->edgeIndex(from, to));
    }

    std::shared_ptr<UnweightedGraph<IsDirected>> topology() const
    {
        return p_topology_;
    }

    size_t vertexNum() const { return p_topology_->vertexNum(); }
    size_t edgeNum() const { return p_topology_->edgeNum(); }
    const std::vector<size_t>& adjacency(size_t v_index) const { return p_topology_->adjacency(v_index); }
    bool validVertex(size_t idx) const { return p_topology_->validVertex(idx); }
    static size_t nullVertex() { return UnweightedGraph<IsDirected>::nullVertex(); }
    const std::map<std::array<size_t, 2>, size_t >& edgeIndices() const { return p_topology_->edgeIndices(); }

protected:
    std::shared_ptr<UnweightedGraph<IsDirected>> p_topology_ = nullptr;
    EdgeProperty<DType> edge_property_;
};


template<typename DType>
using WeightedDirectedGraph = WeightedGraph<DType, true>;
template<typename DType>
using WeightedUnirectedGraph = WeightedGraph<DType, false>;
#endif
} // namespace mxm


#endif // _GRAPH_BASE_H_
