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

template <size_t N>
const Block& rotBlk() { return se::rotBlk<N>(); }

template <size_t N>
const Block& traBlk() { return se::traBlk<N>(); }

template<typename DType>
bool isValid(const Matrix<DType>& mat, DType* error=nullptr, DType tol=eps<DType>())
{
    return true;
}

template <size_t N, typename DType>
Matrix<DType> inv(const Matrix<DType>& mat)
{
    // return mat.T();
    auto ret = Matrix<DType>::identity(N+1);
    auto rot_inv = mat(rotBlk<N>()).T();
    ret(rotBlk<N>()) = rot_inv;
    ret(traBlk<N>()) = -rot_inv.matmul(mat(traBlk<N>()));
    return ret;
}

template <size_t N, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
log(const Matrix<DType>& grp)
{
    auto alg_rt = Matrix<DType>::identity(N+1);

    auto skew = SO::log<N>(grp(rotBlk<N>()));
    alg_rt(rotBlk<N>()) = skew;
    alg_rt(traBlk<N>()) = so::jacobInv(skew).matmul(grp(traBlk<N>()));

    return alg_rt;
}

// denote r = exp(z)
// z = log(r)
// r^y = (e^z)^y = exp(y*z) = exp(y*log(r))
template<size_t N, typename DType>
std::enable_if_t<3 == N,Matrix<DType>>
pow(const Matrix<DType>& pose, DType y)
{
    return se::exp<N>(y * SE::log<N>(pose));
}

template<typename DType>
Matrix<DType> interp(const Matrix<DType>& p0, const Matrix<DType>& p1, DType t)
{
    static const size_t N=3;
    return SE::pow<N>(p1.matmul(SE::inv<N>(p0)), t).matmul(p0);
}

// left lie derivative:
// d Tp / dR
template<size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
derivPoint(const Matrix<DType>& mat, const Matrix<DType>& p_homo)
{
    auto p = p_homo(Block({0, N},{}));
    Matrix<DType> ret = Matrix<DType>::zeros({N + SO::dof<N>(), N + SO::dof<N>()});
    ret(Block({0, N}, {N, end()})) = -so::wedge(mat(rotBlk<N>()).matmul(p) + mat(traBlk<N>()));
    for(size_t i = 0; i < N; i++) ret(i,i) = 1;
    return ret;
}

template<size_t N=3, typename DType>
Matrix<DType> adj(const Matrix<DType> mat)
{
    Matrix<DType> ret = Matrix<DType>::zeros({2*N, 2*N});
    auto rot = mat(rotBlk<N>());
    ret(Block({0, N}, {0, N})) = rot;
    ret(Block({0, N}, {N, 2*N})) = so::wedge(mat(traBlk<N>())).matmul(rot);
    ret(Block({N, 2*N}, {N, 2*N})) = rot;
    return ret;
}


} // namespace SE

// lie algebra of the adjoint of special euclidean algebra
namespace adse
{

} // namespace adse

// lie group of the adjoint of special euclidean group
namespace AdSE
{
template<size_t N=3, typename DType>
Matrix<DType> inv(const Matrix<DType> mat)
{
    Matrix<DType> ret = Matrix<DType>::zeros({2*N, 2*N});
    assert(false);
    return ret;
}
} // namespace AdSE



} // namespace mxm


#endif //_LIE_SPECIAL_EUCLIDEAN_H_