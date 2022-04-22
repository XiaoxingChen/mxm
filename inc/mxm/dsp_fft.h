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
    DType angle_step = -2 * M_PI / static_cast<DType>(n);
    for(size_t i = 0; i < n; i++)
    {
        for(size_t j = i; j < n; j++)
        {
            Complex<DType> angle{0, angle_step*i*j};
            Complex<DType> val = mxm::exp(angle);
            mat(i,j) = val;
            mat(j,i) = val;
        }
    }

    return scale * mat;
}

} // namespace mxm

#endif // __DSP_FFT_H__
