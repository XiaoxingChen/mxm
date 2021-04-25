#if !defined(_GRAPH_BASE_H_)
#define _GRAPH_BASE_H_

#include "linalg_mat.h"
#include <map>

namespace mxm
{

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
} // namespace mxm


#endif // _GRAPH_BASE_H_
