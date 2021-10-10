#if !defined(__DSP_FFT_H__)
#define __DSP_FFT_H__

#include "mxm/linalg.h"

namespace mxm
{

// https://en.wikipedia.org/wiki/DFT_matrix
template <typename DType=float>
Matrix<Complex<DType>> dftMatrix(size_t n)
{
    Matrix<Complex<DType>> mat = Matrix<Complex<DType>>::ones({n,n});
    DType scale = 1. / sqrt(static_cast<DType>(n));
    DType angle = -2 * M_PI / static_cast<DType>(n);
    auto omega = Complex<DType>{cos(angle), sin(angle)};
    auto inc = Matrix<Complex<DType>>::ones({1,n});
    for(size_t i = 1; i < n; i++) inc(0, i) = inc(0, i - 1) * omega;

    for(size_t i = 1; i < n; i++)
    {
        mat(Row(i)) = mat(Row(i - 1)) * inc;
    }
    return scale * mat;
}

} // namespace mxm

#endif // __DSP_FFT_H__
