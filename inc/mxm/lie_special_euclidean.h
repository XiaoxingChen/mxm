#if !defined(_LIE_SPECIAL_EUCLIDEAN_H_)
#define _LIE_SPECIAL_EUCLIDEAN_H_

#include "linalg.h"
#include "lie_special_orthogonal.h"

namespace mxm
{

namespace se
{

template<typename DType>
bool isValid(const Matrix<DType>& mat, DType* error=nullptr, DType tol=eps<DType>())
{
    return true;
}

template <size_t N>
const Block& rotBlk()
{
    static const Block blk({0,N},{0,N});
    return blk;
}

template <size_t N>
const Block& traBlk()
{
    static const Block blk({0,N}, {N, N+1});
    return blk;
}

template <size_t N, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
exp(const Matrix<DType>& alg)
{
    auto mat_rt = Matrix<DType>::identity(N+1);

    auto skew = alg(rotBlk<N>());
    mat_rt(rotBlk<N>()) = so::exp<N>(skew);
    mat_rt(traBlk<N>()) = so::jacob(skew).matmul(alg( traBlk<N>()));

    return mat_rt;
}

} // namespace se

namespace SE
{

template<typename DType>
bool isValid(const Matrix<DType>& mat, DType* error=nullptr, DType tol=eps<DType>())
{
    return true;
}

template <size_t N, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
log(const Matrix<DType>& grp)
{
    auto alg_rt = Matrix<DType>::identity(N+1);

    auto skew = SO::log<N>(grp(se::rotBlk<N>()));
    alg_rt(se::rotBlk<N>()) = skew;
    alg_rt(se::traBlk<N>()) = so::jacobInv(skew).matmul(grp(se::traBlk<N>()));

    return alg_rt;
}

} // namespace SE

} // namespace mxm


#endif //_LIE_SPECIAL_EUCLIDEAN_H_