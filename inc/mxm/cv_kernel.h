#if !defined(_CV_KERNEL_H_)
#define _CV_KERNEL_H_
#include "linalg_mat.h"

namespace mxm
{

namespace kernel
{

template<typename DType=float>
typename std::enable_if_t<std::is_floating_point<DType>::value, const Matrix<DType>&>
sobel()
{
    const static Matrix<DType> kernel({3,3}, {
        -0.125, -0.25, -0.125,
        0, 0, 0,
        0.125, 0.25, 0.125});
    return kernel;
}

template<typename DType=float>
typename std::enable_if_t<std::is_floating_point<DType>::value, Matrix<DType>>
gauss(size_t size)
{
    if(size == 1) return Matrix<DType>({1,1},{1});
    size_t layer = size - 1;
    Matrix<DType> vec = Matrix<DType>({size, 1});
    for(size_t i = 0; i < (size + 1)/ 2; i++)
    {
        float comb = mxm::combinations(layer, i);
        vec(i,0) = comb;
        vec(size - i - 1,0) = comb;
    }
    // std::cout << "vec: " << vec.T().str() << std::endl;
    vec *= (DType(1) / (1 << layer));

    return vec.matmul(vec.T());
}

} // namespace kernel
} // namespace mxm


#endif // _CV_KERNEL_H_
