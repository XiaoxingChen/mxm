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

template <size_t N=3>
const Block& rotBlk()
{
    static const Block blk({0,N},{0,N});
    return blk;
}

template <size_t N=3>
const Block& traBlk()
{
    static const Block blk({0,N}, {N, N+1});
    return blk;
}

template <size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
exp(const Matrix<DType>& alg)
{
    auto mat_rt = Matrix<DType>::identity(N+1);

    auto skew = alg(rotBlk<N>());
    mat_rt(rotBlk<N>()) = so::exp<N>(skew);
    mat_rt(traBlk<N>()) = so::jacob(skew).matmul(alg( traBlk<N>()));

    return mat_rt;
}

template <size_t N=3, typename DeriveType>
std::enable_if_t<3 == N, Matrix<typename Traits<DeriveType>::EntryType>>
wedge(const MatrixBase<DeriveType>& vec)
{
    using DType = typename Traits<DeriveType>::EntryType;
    auto mat_se = Matrix<DType>::zeros({N+1, N+1});
    mat_se(traBlk<N>()) = vec(Block({0, N},{}));
    mat_se(rotBlk<N>()) = so::wedge(vec(Block({N, end()},{})));
    return mat_se;
}

template <size_t N=3, typename DeriveType>
std::enable_if_t<3 == N, Vector<typename Traits<DeriveType>::EntryType>>
vee(const MatrixBase<DeriveType>& mat_se)
{
    using DType = typename Traits<DeriveType>::EntryType;
    Vector<DType> vec(N*2);
    vec(Block({0, N},{})) = mat_se(traBlk<N>());
    vec(Block({N, end()},{})) = so::vee(mat_se(rotBlk<N>()));
    return vec;
}

// formula 7.51a
template <size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
adj(const Matrix<DType>& alg)
{
    Matrix<DType> ret = Matrix<DType>::identity(2*N);
    auto rot = alg(rotBlk<N>());
    ret.setBlock(0,0, rot);
    ret.setBlock(N,N, rot);
    ret.setBlock(0,N, so::wedge(alg(traBlk<N>())));
    return ret;
}


// formula 7.86a
template <size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
matQ(const Matrix<DType>& mat)
{
    Matrix<DType> phi = mat(rotBlk<N>());
    DType angle = so::findAngle<N>(phi);
    if(abs(angle) < eps<DType>()) return Matrix<DType>::zeros({N,N});

    DType i_angle = DType(1) / angle;
    DType i_angle3 = i_angle * i_angle * i_angle;
    DType i_angle4 = i_angle * i_angle3;
    DType i_angle5 = i_angle * i_angle4;


    Matrix<DType> phi2 = phi.matmul(phi);
    Matrix<DType> rho = so::wedge(mat(traBlk<N>()));
    Matrix<DType> phi_rho_phi = phi.matmul(rho).matmul(phi);

    Matrix<DType> mat_q = rho * 0.5;
    mat_q += (angle - sin(angle)) * i_angle3 * (phi.matmul(rho) + rho.matmul(phi) + phi_rho_phi);
    mat_q += (angle * angle + 2 * cos(angle) - 2) * 0.5 * i_angle4 * (phi2.matmul(rho) + rho.matmul(phi2) - phi_rho_phi * 3.);
    mat_q += (2 * angle - 3 * sin(angle) + angle * cos(angle)) * 0.5 * i_angle5 * (phi_rho_phi.matmul(phi) + phi.matmul(phi_rho_phi));
    return mat_q;
}

template <size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
jacob(const Matrix<DType>& mat)
{
    auto ret = Matrix<DType>::identity(2*N);
    auto so_jac = so::jacob(mat(rotBlk<N>()));
    ret.setBlock(0,0, so_jac);
    ret.setBlock(N,N, so_jac);
    ret.setBlock(0,N, matQ(mat));
    return ret;
}

// formula 7.95b
template <size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
jacobInv(const Matrix<DType>& mat)
{
    auto ret = Matrix<DType>::identity(2*N);
    auto so_jac_inv = so::jacobInv(mat(rotBlk<N>()));
    ret.setBlock(0,0, so_jac_inv);
    ret.setBlock(N,N, so_jac_inv);
    ret.setBlock(0,N, - so_jac_inv.matmul(matQ(mat)).matmul(so_jac_inv) );
    return ret;
}

// formula 7.94
template <size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
jacob2(const Matrix<DType>& mat)
{

}

} // namespace se

namespace SE
{

template <size_t N=3>
const Block& rotBlk() { return se::rotBlk<N>(); }

template <size_t N=3>
const Block& traBlk() { return se::traBlk<N>(); }

template<typename DType>
bool isValid(const Matrix<DType>& mat, DType* error=nullptr, DType tol=eps<DType>())
{
    return true;
}

template <size_t N=3, typename DType>
Matrix<DType> inv(const Matrix<DType>& mat)
{
    // return mat.T();
    auto ret = Matrix<DType>::identity(N+1);
    auto rot_inv = mat(rotBlk<N>()).T();
    ret(rotBlk<N>()) = rot_inv;
    ret(traBlk<N>()) = -rot_inv.matmul(mat(traBlk<N>()));
    return ret;
}

template <size_t N=3, typename DType>
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
    ret.setBlock(0,0, rot);
    ret.setBlock(N,N, rot);
    ret.setBlock(0,N, so::wedge(mat(traBlk<N>())).matmul(rot));

    return ret;
}

// left lie derivative:
// d ln(T T_1) / dT
// formulate 7.83a
template<size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
derivDistance(const Matrix<DType>& tf, const Matrix<DType>& tf_1)
{
    return se::jacobInv<N>(SE::log<N>(tf.matmul(tf_1)));
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