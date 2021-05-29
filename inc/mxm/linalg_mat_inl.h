#if !defined(_LINALG_MAT_INL_H_)
#define _LINALG_MAT_INL_H_

#if defined(MXM_COMPILED_LIB)
#include "linalg_mat.h"
#endif // MXM_COMPILED_LIB


namespace mxm
{

AutoShape fixRow(size_t n) { return AutoShape(n, 1, true, false); }
AutoShape fixCol(size_t n) { return AutoShape(0, n, false, true); }

std::array<size_t, 2> AutoShape::deduct(size_t total_num) const
{
    if(defined_[0] && defined_[1]) return shape_;
    std::array<size_t, 2> ret(shape_);
    if(!defined_[0]) ret[0] = total_num / shape_[1];
    if(!defined_[1]) ret[1] = total_num / shape_[0];
    return ret;
}

void updateOffset(
    Shape& abs_offset, const Shape& inc_offset, bool same_major)
{
    if(same_major)
    {
        abs_offset[0] += inc_offset[0];
        abs_offset[1] += inc_offset[1];
    }else
    {
        abs_offset[0] += inc_offset[1];
        abs_offset[1] += inc_offset[0];
    }
}

} // namespace mxm

#endif // _LINALG_MAT_INL_H_
