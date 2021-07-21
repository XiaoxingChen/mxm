#if !defined(_LINALG_MAT_INL_H_)
#define _LINALG_MAT_INL_H_

#if defined(MXM_COMPILED_LIB)
#include "linalg_mat.h"
#endif // MXM_COMPILED_LIB


namespace mxm
{
#if 0
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
#else
AutoShape fixRow(size_t n) { return AutoShape(AutoShape::eRowDefined, n, 0); }
AutoShape fixCol(size_t n) { return AutoShape(AutoShape::eColDefined, 0, n); }
AutoShape square() { return AutoShape(AutoShape::eSquare, 0, 0); }

std::array<size_t, 2> AutoShape::deduct(size_t total_num) const
{
    if(eFullyDefined == state_) return shape_;
    std::array<size_t, 2> ret(shape_);
    if(eColDefined == state_) ret[0] = total_num / shape_[1];
    else if(eRowDefined == state_) ret[1] = total_num / shape_[0];
    else if(eSquare == state_)
    {
        size_t n(static_cast<size_t>(sqrt(total_num) + 0.5));
        ret = {n,n};
    }
    return ret;
}
#endif

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
