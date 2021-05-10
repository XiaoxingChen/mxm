#if !defined(_LINALG_MAT_INL_H_)
#define _LINALG_MAT_INL_H_

#if defined(MXM_COMPILED_LIB)
#include "linalg_mat.h"
#endif // MXM_COMPILED_LIB


namespace mxm
{

AutoShape fixRow(size_t n) { return AutoShape(n, -1); }
AutoShape fixCol(size_t n) { return AutoShape(-1, n); }

std::array<size_t, 2> AutoShape::deduct(size_t total_num) const
{
    std::array<size_t, 2> ret;
    if(raw_shape_[0] < 0)
    {
        ret[0] = total_num / raw_shape_[1];
        ret[1] = raw_shape_[1];
    }else if(raw_shape_[1] < 0)
    {
        ret[0] = raw_shape_[0];
        ret[1] = total_num / raw_shape_[0];
    }else
    {
        ret[0] = raw_shape_[0];
        ret[1] = raw_shape_[1];
    }

    if(ret[0] * ret[1] != total_num && total_num != 0)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    return ret;

}

} // namespace mxm

#endif // _LINALG_MAT_INL_H_
