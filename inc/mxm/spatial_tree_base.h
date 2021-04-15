#if !defined(__SPATIAL_TREE_BASE_H__)
#define __SPATIAL_TREE_BASE_H__
#include <array>
namespace mxm
{

struct RangeNode
{
    size_t node_idx;
    std::array<size_t, 2> range;
};

} // namespace mxm


#endif // __SPATIAL_TREE_BASE_H__
