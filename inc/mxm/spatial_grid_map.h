#if !defined(__SPATIAL_GRID_MAP_H__)
#define __SPATIAL_GRID_MAP_H__

#include "linalg.h"

namespace mxm
{

template<typename DType>
Matrix<uint32_t> quantize(const Matrix<DType>& pts, const Matrix<DType>& offset, DType scale)
{
    Matrix<uint32_t> result(pts.shape());
    result.traverse([&](auto i, auto j){
        result(i,j) = static_cast<uint32_t> ((pts(i,j)  + offset(i,0)) * scale + 0.5);
    });
    return result;
}

template<typename DType>
Matrix<uint32_t> quantize(const Matrix<DType>& pts, uint32_t resolution)
{
    Matrix<DType> bounds = boundary(pts);
    Vector<DType> shape = bounds(Col(1)) - bounds(Col(0));
    DType scale = static_cast<DType>(resolution - 1) / mxm::max(shape);
    return quantize(pts, bounds(Col(0)) * -1, scale);
}


} // namespace mxm


#endif // __SPATIAL_GRID_MAP_H__
