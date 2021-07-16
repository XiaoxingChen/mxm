#if !defined(_LINALG_SOLVE_H)
#define _LINALG_SOLVE_H

#include "mxm/linalg_mat_ref.h"
#include "mxm/linalg_complex.h"
#include "mxm/linalg_norm.h"
#include "mxm/string.h"
#include <iostream>
#include <algorithm>


namespace mxm
{

template<typename DeriveType1, typename DeriveType2>
Matrix<typename Traits<DeriveType1>::EntryType>
solveLUTriangle(const MatrixBase<DeriveType1>& mat_in, const MatrixBase<DeriveType2>& b_in, bool l_tri)
{
    using DType = typename Traits<DeriveType1>::EntryType;
    auto & mat = reinterpret_cast<const DeriveType1&>(mat_in);
    auto & b = reinterpret_cast<const DeriveType2&>(b_in);

    if(!mat.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    if(mat.shape(0) != b.shape(0))
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    size_t idx_start = l_tri ? 0 : b.shape(0) - 1;
    int step = l_tri ? 1 : - 1;
    Matrix<DType> x(Matrix<DType>::zeros(b.shape()));

    for(size_t w = 0; w < x.shape(1); w++)
    {
        for(size_t i = idx_start, _1 = 0; _1 < b.shape(0); _1 ++, i += step)
        {
            DType off_diag_sum(0);
            for(size_t j = idx_start, _2 = 0; _2 < _1; _2 ++, j += step)
            {
                off_diag_sum += mat(i, j) * x(j,w);
            }
            x(i,w) = (b(i,w) - off_diag_sum) * inv(mat(i,i));
        }
    }

    return x;
}


template<typename DType>
Matrix<DType> solveLowerTriangle(const Matrix<DType>& lower, const Matrix<DType>& b)
{
    return solveLUTriangle(lower, b, 1);
}

template<typename DeriveType1, typename DeriveType2>
Matrix<typename Traits<DeriveType1>::EntryType>
solveUpperTriangle(const MatrixBase<DeriveType1>& upper, const MatrixBase<DeriveType2>& b)
{
    return solveLUTriangle(upper, b, 0);
}

namespace qr
{
// u: column vector, direction
// a: column vectors
template<typename DeriveType1, typename DeriveType2>
typename Traits<DeriveType1>::DerefType
project(const MatrixBase<DeriveType1>& u, const MatrixBase<DeriveType2>& a)
{
    auto u_dir = u.normalized();
    return u_dir.matmul(u_dir.T()).matmul(a);
}

inline std::vector<std::array<size_t, 2>> upperTrianglizeSequence(size_t cols)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < cols; j++)
        for(size_t i = cols-1; i > j; i--)
            ret.push_back({i,j});
    return ret;
}
inline std::vector<std::array<size_t, 2>> upperHessenbergizeSequence(size_t cols)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < cols; j++)
        for(size_t i = cols-1; i > j+1; i--)
            ret.push_back({i,j});
    return ret;
}
inline std::vector<std::array<size_t, 2>> subdiagonalSequence(const Shape& shape)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < shape[1]-1; j++)
        ret.push_back({j+1,j});
    if(shape[0] > shape[1]) ret.push_back({shape[1], shape[1]-1});
    return ret;
}
inline std::vector<std::array<size_t, 2>> superdiagonalSequence(const Shape& shape)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 1; j < shape[1]; j++)
        ret.push_back({j-1,j});
    return ret;
}
enum TraverseSeq{
    eUpperTrianglize,
    eUpperHessenbergize,
    eSubdiagonal
};
#if 0
template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, Matrix<DType>>::type
givensRotation(const Matrix<DType>& v2)
{
    auto norm_v = v2.norm();
    Matrix<DType> so2({2,2},{
        v2(0,0), v2(1,0),
        -v2(1,0), v2(0,0) });
    so2 *= (decltype(norm_v)(1) / norm_v);
    return so2;
}


template<typename DType, unsigned int N>
Matrix<Hypercomplex<DType, N>>
givensRotation(const Matrix<Hypercomplex<DType, N>>& v2)
{
    auto norm_v = v2.norm();
    Matrix<Hypercomplex<DType, N>> su2({2,2},{
        v2(0,0).conj(), v2(1,0).conj(),
        -v2(1,0), v2(0,0) });
    su2 *= (decltype(norm_v)(1) / norm_v);
    return su2;
}
#else

template<template <class> class MatrixType, class EntryType>
std::enable_if_t<std::is_floating_point<EntryType>::value , Matrix<EntryType>>
givensRotation(const MatrixBase<MatrixType<EntryType>>& v2)
{
    auto norm_v = v2.norm();
    Matrix<EntryType> so2({2,2},{
        v2(0,0), v2(1,0),
        -v2(1,0), v2(0,0) });
    so2 *= (decltype(norm_v)(1) / norm_v);
    return so2;
}

template<template <class> class MatrixType, class EntryType>
std::enable_if_t<IsHypercomplexMatrix<MatrixType<EntryType>>::value , Matrix<EntryType>>
givensRotation(const MatrixBase<MatrixType<EntryType>>& v2)
{
    auto norm_v = v2.norm();
    Matrix<EntryType> su2({2,2},{
        v2(0,0).conj(), v2(1,0).conj(),
        -v2(1,0), v2(0,0) });
    su2 *= (decltype(norm_v)(1) / norm_v);
    return su2;
}

#endif

//
// idx_seq: TraverseSeq
// symmetric: output will be QHQ' = mat if enabled.
// output: {Q, R}
// Reference:
// [1] https://www.math.usm.edu/lambers/mat610/sum10/lecture9.pdf
template<typename DeriveType>
std::array<Matrix<typename Traits<DeriveType>::EntryType>, 2>
decomposeByRotation(const MatrixBase<DeriveType>& mat_in, TraverseSeq idx_seq=eUpperTrianglize, bool symmetric=false)
{
    using DType = typename Traits<DeriveType>::EntryType;
    if(!mat_in.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    std::array<Matrix<DType>, 2> ret;
    Matrix<DType>& mat(ret[1]);
    mat = mat_in;
    const size_t& n = mat.shape(0);
    Matrix<DType>& rot(ret[0]);
    rot = Matrix<DType>::identity(n);

    std::vector<std::array<size_t, 2>> seq;
    if(idx_seq == eUpperTrianglize) seq = upperTrianglizeSequence(mat_in.shape(1));
    else if (idx_seq == eUpperHessenbergize) seq = upperHessenbergizeSequence(mat_in.shape(1));
    else if (idx_seq == eSubdiagonal)  seq = subdiagonalSequence(mat_in.shape());

    for(auto& idx: seq)
    {
        auto i = idx[0];
        auto j = idx[1];
        if(norm(mat(i,j)) < eps() * eps()) { continue; }

        Matrix<DType> sub_rot = Matrix<DType>::identity(n);

        auto so2 = givensRotation(mat(Block({i-1, i+1}, {j, j+1})));

        sub_rot.setBlock(i-1, i-1, so2);

        rot = sub_rot.matmul(rot);
        mat = sub_rot.matmul(mat);
        if(symmetric)
            mat = mat.matmul(sub_rot.T());

        // if(std::is_same<DType, Complex<FloatType>>::value)
        // {
        //     // std::cout << "mat: \n" << mxm::to_string(mat) << std::endl;
        //     std::cout << "recover: \n" << mxm::to_string((conj(rot.T())).matmul(mat)) << std::endl;
        // }
    }
    rot = conj(rot.T());
    return ret;
}

template<typename DType>
Matrix<DType> calcMatQFromReflection(const Matrix<DType>& mat)
{
    //reference: https://en.wikipedia.org/wiki/QR_decomposition
    if(!mat.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    Matrix<DType> mat_u({mat.shape(0), mat.shape(1)});
    for(size_t i = 0; i < mat.shape(0); i++)
    {
        mat_u(Col(i)) = mat(Col(i));
        for(size_t j = 0; j < i; j++)
        {
            mat_u(Col(i)) -= project(mat_u(Col(j)), mat(Col(i)));
        }
    }

    // std::cout << "Bq: \n" << mxm::to_string(mat_u);
    for(size_t i = 0; i < mat.shape(0); i++)
    {
        mat_u(Col(i)) = mat_u(Col(i)).normalized();
    }

    return mat_u;
}

template<typename DeriveType>
auto calcMatQ(const MatrixBase<DeriveType>& mat)
{
    return decomposeByRotation(mat)[0];
}

template<typename DType>
Matrix<DType> solve(const Matrix<DType>& mat_a, const Matrix<DType>& b)
{
    // Matrix<DType> mat_q(calcMatQ(mat_a));
    if( norm(mat_a) < eps())
    {
        std::cout << "Zero Matrix!" << std::endl;
        return Matrix<DType>::zeros(b.shape());
    }
    auto q_r = decomposeByRotation(mat_a);

    if((q_r[0].matmul(conj(q_r[0].T())) - Matrix<DType>::identity(q_r[0].shape(0)) ).norm()  > 10 * eps())
    {
        // singular
        auto tmp = q_r[0].matmul(conj(q_r[0].T()));
        std::cout << "q.matmul(q.T()): \n" << mxm::to_string(tmp) << std::endl;
        std::cout << "norm: " << (tmp - Matrix<DType>::identity(q_r[0].shape(0)) ).norm() << std::endl;
        std::cout << "Singular Matrix!" << std::endl;
        return Matrix<DType>::zeros(b.shape());
    }

    // Matrix<DType> mat_r(mat_q.T().matmul(mat_a));
    Matrix<DType> x(solveUpperTriangle(q_r[1], conj(q_r[0].T()).matmul(b)));
    return x;
}

// Reference:
// https://math.stackexchange.com/questions/1262363/convergence-of-qr-algorithm-to-upper-triangular-matrix
template<typename DType>
DType wilkinsonShiftStrategy(const Matrix<DType> mat_a)
{
    size_t n = mat_a.shape(1);
    DType sigma = 0.5 * (mat_a(n-2,n-2) - mat_a(n-1, n-1));
    if(abs(sigma) < std::numeric_limits<DType>::epsilon() && abs(mat_a(n-1, n-1)) < std::numeric_limits<DType>::epsilon())
        return DType(0);
    DType sign = sigma > 0 ? 1. : -1.;
    DType mu = mat_a(n-1, n-1) - (sign * mat_a(n-1, n-2) * mat_a(n-1, n-2)) / (abs(sigma) + sqrt(sigma * sigma) + mat_a(n-1, n-1) * mat_a(n-1, n-1));
    return mu;
}

template<typename DType>
DType errorOrthogonalBlockDiagonal(const Matrix<DType> mat_a)
{
    DType sum = 0.;
    auto seq = superdiagonalSequence(mat_a.shape());
    for(auto & idx : seq)
    {
        DType res = mat_a(idx[0], idx[1]) + mat_a(idx[1], idx[0]);
        sum += res * res;
    }
    return sum;
}

} // namespace qr
#if 0
template <typename DType>
DType Matrix<DType>::det() const
{
    if(square() && shape(0) == 2)
        return (*this)(0,0) * (*this)(1,1) - (*this)(1,0)*(*this)(0,1);

    Matrix<DType> mat_q(qr::calcMatQ(*this));

    // check full rank
    if((mat_q.matmul(mat_q)).trace() - shape(0) > eps())
        return 0.;

    Matrix<DType> mat_r(mat_q.T().matmul(*this));
    DType det(1);
    for(size_t i = 0; i < shape(0); i++) det *= mat_r(i,i);
    return det;
}
#endif
template <typename DeriveType>
typename Traits<DeriveType>::EntryType
det(const MatrixBase<DeriveType>& mat_in)
{
    using DType = typename Traits<DeriveType>::EntryType;
    auto & mat = reinterpret_cast<const DeriveType&>(mat_in);
    if(mat.square() && mat.shape(0) == 2)
        return mat(0,0) * mat(1,1) - mat(1,0)*mat(0,1);

    Matrix<DType> mat_q(qr::calcMatQ(mat));

    // check full rank
    if((mat_q.matmul(mat_q)).trace() - mat.shape(0) > eps())
        return 0.;

    Matrix<DType> mat_r(mat_q.T().matmul(mat));
    DType det(1);
    for(size_t i = 0; i < mat.shape(0); i++) det *= mat_r(i,i);
    return det;
}
#if 0
template <typename DType>
Matrix<DType> Matrix<DType>::inv() const
{
    if(!square()) return Matrix<DType>::zeros(shape());
    return qr::solve(*this, Matrix<DType>::identity(shape(0)));
}
#endif

// Inversion
template<typename DeriveType, typename>
Matrix<typename Traits<DeriveType>::EntryType>
inv(const MatrixBase<DeriveType>& mat)
{
    using EntryType = typename Traits<DeriveType>::EntryType;
    using ArithType = typename Traits<DeriveType>::ArithType;
    auto& self = reinterpret_cast<const DeriveType&>(mat);
    if(!mat.square()) return Matrix<EntryType>::zeros(self.shape());
    return qr::solve(self, Matrix<EntryType>::identity(self.shape(0)));
}

template <template <class> class MatrixType, typename DType>
std::vector<Complex<DType>> eigvals2x2(const MatrixBase<MatrixType<DType>>& mat)
{
    std::vector<Complex<DType>> ret(2);
    DType tr = mat.trace();
    DType det = mxm::det(mat);
    DType delta = tr*tr - 4 * det;
    if(delta >= 0)
    {
        ret.at(0) = Complex<DType>({DType(0.5) * (tr + sqrt(delta)), 0});
        ret.at(1) = Complex<DType>({DType(0.5) * (tr - sqrt(delta)), 0});
        // std::cout << "tr: " << tr <<", det: " << det << ", delta: " << delta << std::endl;
        // std::cout << "mat: \n" << mxm::to_string(mat) << std::endl;
    }else
    {
        ret.at(0) = Complex<DType>({DType(0.5) * tr,  DType(0.5) * sqrt(-delta)});
        ret.at(1) = Complex<DType>({DType(0.5) * tr, -DType(0.5) * sqrt(-delta)});
        // std::cout << "<0" << std::endl;
    }
    return ret;
}

template<typename DType>
Matrix<DType> shiftedQRIteration(
    const Matrix<DType>& mat,
    Matrix<DType>* p_orthogonal=nullptr,
    qr::TraverseSeq idx_seq=qr::eUpperTrianglize,
    size_t max_it=40,
    DType tol=std::numeric_limits<DType>::epsilon() * std::numeric_limits<DType>::epsilon())
{
    std::array<Matrix<DType>, 2> q_r;
    Matrix<DType> ret(mat);
    Matrix<DType> eye = Matrix<DType>::identity(mat.shape(0));

    for(size_t i = 0; i < max_it; i++)
    {
        // FloatType rho = 1;
        DType rho = qr::wilkinsonShiftStrategy(ret);
        Matrix<DType> shift = eye * rho;
        q_r = qr::decomposeByRotation(ret - shift, idx_seq);
        ret = q_r[1].matmul(q_r[0]) + shift;
        if(p_orthogonal) *p_orthogonal = p_orthogonal->matmul(q_r[0]);
        if(qr::errorOrthogonalBlockDiagonal(q_r[0]) < tol) break;
    }
    return ret;
}

// A = QTQ', where
// Q is a real, orthogonal matrix,
// T is a real, quasi-upper-triangular matrix that has a block upper-triangular structure
// Reference:
// [1] http://www.math.usm.edu/lambers/mat610/sum10/lecture15.pdf
template<typename DType>
Matrix<DType> realSchurDecomposition(const Matrix<DType>& mat, Matrix<DType>* p_orthogonal=nullptr)
{
    size_t n = mat.shape(0);
    if(n < 3)
    {
        if(p_orthogonal) *p_orthogonal = Matrix<DType>::identity(n);
        return mat;
    }
#if 0
    const size_t max_it = 30;
    const DType tol = std::numeric_limits<DType>::epsilon() * std::numeric_limits<DType>::epsilon();

    auto q_hesson = qr::decomposeByRotation(mat, qr::eUpperHessenbergize, true);
    Matrix<DType> block_diag = q_hesson[1];
    if(p_orthogonal) *p_orthogonal = q_hesson[0];

    std::array<Matrix<DType>, 2> q_r;

    for(size_t i = 0; i < max_it; i++)
    {
        // FloatType rho = 1;
        DType rho = qr::wilkinsonShiftStrategy(block_diag);
        Matrix<DType> shift = Matrix<DType>::identity(n) * rho;
        q_r = qr::decomposeByRotation(block_diag - shift, qr::eSubdiagonal);
        block_diag = q_r[1].matmul(q_r[0]) + shift;
        if(p_orthogonal) *p_orthogonal = p_orthogonal->matmul(q_r[0]);
        if(qr::errorOrthogonalBlockDiagonal(q_r[0]) < tol) break;
    }

    return block_diag;
#else
    auto q_hesson = qr::decomposeByRotation(mat, qr::eUpperHessenbergize, true);
    if(p_orthogonal) *p_orthogonal = q_hesson[0];
    return shiftedQRIteration(q_hesson[1], p_orthogonal, qr::eSubdiagonal);
#endif
}

// Reference:
// http://www.math.usm.edu/lambers/mat610/sum10/lecture15.pdf
// https://www.cs.purdue.edu/homes/skeel/CS515/4.pdf
template <typename DType>
std::vector<Complex<DType>> eigvals(const Matrix<DType>& mat)
{
    if(!mat.square()) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    size_t n = mat.shape(0);
    auto quasi_upper_triangle = realSchurDecomposition(mat);
    // std::cout << mxm::to_string(quasi_upper_triangle, 12) << std::endl;
    std::vector<Complex<DType>> ret(n);
    for(size_t i = 0; i < n;)
    {
        // if(i < n-1 && abs(q_r[0](i,i) - q_r[0](i+1, i+1)) < tol && abs(q_r[0](i+1,i)) > tol)
        if(i < n-1 && abs(quasi_upper_triangle(i+1,i)) > 50 * std::numeric_limits<DType>::epsilon())
        {
            // ret.setBlock(i,0, eigvals2x2(quasi(Block({i,i+2},{i,i+2}))));
            auto eig_pair = eigvals2x2(quasi_upper_triangle(Block({i,i+2},{i,i+2})));
            ret.at(i) = eig_pair.at(0);
            ret.at(i+1) = eig_pair.at(1);
            i += 2;
        }else
        {
            ret.at(i) = Complex<DType>({quasi_upper_triangle(i,i), 0});
            i++;
        }
    }

    return ret;
}

// DType is required to be double or float
// std::is_float_point<DType>::value == true
// todo: consider reducing `max_it` according to [2]
// Reference:
// [1] https://dspace.mit.edu/bitstream/handle/1721.1/75282/18-335j-fall-2006/contents/lecture-notes/lec16handout6pp.pdf
// [2] https://www5.in.tum.de/lehre/vorlesungen/konkr_math/WS_11_12/vorl/NumPro_WS1112_Vorlesung_Kapitel_7.pdf#page=18
template <typename DType>
std::array<Matrix<DType>, 2> symmetricEig(const Matrix<DType>& mat)
{
    if(!mat.square()) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    size_t n = mat.shape(0);
    DType tol = 1e-13;
    size_t max_it = 100;

    std::array<Matrix<DType>, 2> val_vec{Matrix<DType>({n,1}), Matrix<DType>({n,n})};

    auto eig_vecs = Matrix<DType>::identity(n);

    Matrix<DType> mat_k = mat;
    if(n > 2)
    {
        auto q_r = qr::decomposeByRotation(mat, qr::eUpperHessenbergize, true);
        eig_vecs = q_r[0].matmul(eig_vecs);
        mat_k = q_r[1];
    }

    // std::cout << mxm::to_string(mat_k) << std::endl;

    for(size_t i = 0; i < max_it; i++)
    {
        DType rho = qr::wilkinsonShiftStrategy(mat_k);
        Matrix<DType> shift = Matrix<DType>::identity(n) * rho;
        auto q_r = qr::decomposeByRotation(mat_k - shift, qr::eSubdiagonal);
        mat_k = q_r[1].matmul(q_r[0]) + shift;
        eig_vecs = eig_vecs.matmul(q_r[0]);
        if(qr::errorOrthogonalBlockDiagonal(mat_k) < tol) break;
        // std::cout << "i: " << i << ", err: " << qr::errorOrthogonalBlockDiagonal(mat_k) << std::endl;
    }

    // std::cout << mxm::to_string(mat_k) << std::endl;
    std::vector<size_t> idx_buffer;
    for(size_t i = 0; i < n; i++) idx_buffer.push_back(i);
    std::sort(idx_buffer.begin(), idx_buffer.end(), [&](size_t i, size_t j) { return mat_k(i,i) > mat_k(j,j); });

    // std::cout << mxm::to_string(idx_buffer) << std::endl;

    for(size_t i = 0; i < n; i++)
    {
        size_t idx = idx_buffer.at(i);
        val_vec[0](i,0) = mat_k(idx,idx);
        val_vec[1](Col(i)) = eig_vecs(Col(idx));
    }

    return val_vec;
}

// References:
// https://en.wikipedia.org/wiki/Inverse_iteration
template<typename DType>
Matrix<DType> inverseIteration(
    const Matrix<DType>& mat_a,
    const DType& eigenvalue,
    const Matrix<DType>& guess)
{
    Matrix<DType> bk = guess;
    size_t n = guess.shape(0);
    typename Traits<DType>::ArithType tol = eps() * 20;
    size_t max_it = 20;
    for(size_t i = 0; i < max_it; i++)
    {
        bk = qr::solve(mat_a - eigenvalue * Matrix<DType>::identity(n), bk);
        bk.normalize();
        auto err = (mat_a.matmul(bk) - eigenvalue * bk).norm();
        // std::cout << "bk: " << mxm::to_string(bk.T()) << ", err: " << err << std::endl;
        if(err < tol) break;
    }
    return bk;
}

template<typename DType>
std::array<Matrix<Complex<DType>>, 2>
eig(const Matrix<DType>& mat)
{
    if(!mat.square()) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    size_t n = mat.shape(0);
    std::array<Matrix<Complex<DType>>, 2> ret;
    auto & eig_vals = ret[0];
    auto & eig_vecs = ret[1];
    eig_vals = Matrix<Complex<DType>>({n,1});
    eig_vecs = Matrix<Complex<DType>>({n,n});

    Matrix<Complex<DType>> cmat(mat);

    std::vector<Complex<DType>> prio_que = eigvals(mat);
    std::sort(prio_que.begin(), prio_que.end(), [](auto lhs, auto rhs){return lhs.norm() > rhs.norm();});
    Matrix<Complex<DType>> guess = Matrix<Complex<DType>>::zeros({n,1});
    guess(0,0) = 1;
    for(size_t i = 0;i < prio_que.size(); i++)
    {
        eig_vals(i,0) = prio_que.at(i);
        eig_vecs.setBlock(0,i, inverseIteration(cmat, prio_que.at(i), guess));
    }
    // std::cout << mxm::to_string(ret) << std::endl;
    return ret;
}

#if 1
template<typename DType>
std::array<Matrix<DType>, 3>
svd(const Matrix<DType>& mat)
{
    std::array<Matrix<DType>, 3> u_s_vh;

    size_t n = std::min(mat.shape(0), mat.shape(1));
    u_s_vh[1] = Matrix<DType>({n,1});
    auto inv_sv = Matrix<DType>({n,1});

    std::vector<size_t> idx_buffer(n);
    for(size_t i = 0; i < n; i++) idx_buffer.push_back(i);

    {
        auto val_vec = symmetricEig(mat.matmul(mat.T()));
        u_s_vh[0] = val_vec[1];
        for(size_t i = 0; i < n; i++)
        {
            u_s_vh[1](i,0) = sqrt(val_vec[0](i,0));
            inv_sv(i,0) = DType(1.)/u_s_vh[1](i,0);
        }

        u_s_vh[2] = diagonalMatrix(inv_sv).matmul(u_s_vh[0].T()).matmul(mat);
        // std::cout << mxm::to_string(val_vec[0]) << std::endl;
    }

    {
        // auto val_vec = symmetricEig(mat.T().matmul(mat));
        // u_s_vh[2] = val_vec[1].T();
        // std::cout << mxm::to_string(val_vec[0]) << std::endl;
    }
    return u_s_vh;
}
#endif

// References:
// https://www.cs.cornell.edu/courses/cs6210/2010fa/A6/A6.pdf
template <typename DType>
std::array<Matrix<DType>, 2> tridiagonalizeSkewSymmetric(const Matrix<DType>& skew)
{
    assert(skew.square() && "Skew symmetric matrix must be square!");
#if 0 // by reflection
    const size_t& n = skew.shape(0);
    mat_orth = Matrix<DType>::identity(n);
    for(size_t i = 0; i < n-1; i++)
    {
        Matrix<DType> v = skew(Block({i+1, end()}, {i, i+1}));
        v(0,0) += v.norm();
        v.normalize();
        Matrix<DType> mat_p = Matrix<DType>::identity(n);
        // std::cout << mat_p(Block({i+1, end()}, {i+1, end()})).shape(0) << "," << v.shape(0) <<  std::endl;
        mat_p(Block({i+1, end()}, {i+1, end()})) -= 2* v.matmul(v.T());
        mat_orth = mat_p.matmul(mat_orth);
        skew = mat_p.matmul(skew).matmul(mat_p);
    }
#else // by rotation

    auto ret = qr::decomposeByRotation(skew, qr::eUpperHessenbergize, /*symmetric*/true);
    return ret;
#endif
}

template <typename DType>
std::array<Matrix<DType>, 2> blockDiagonalizeSkewSymmetric(const Matrix<DType>& skew)
{
    if(!skew.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    std::array<Matrix<DType>, 2> ret;
    Matrix<DType>& mat(ret[1]);
    mat = skew;
    const size_t& n = mat.shape(0);
    Matrix<DType>& rot(ret[0]);
    rot = Matrix<DType>::identity(n);

    for(size_t j = 0; j < n; j++)
    {
        for(size_t i = n-1; i > (0 == j % 2 ? j+1 : j); i--)
        {
            std::cout << "(i,j): " << i << "," << j << std::endl;
            DType theta = -atan2(mat(i,j), mat(i-1,j));
            // std::cout << "theta: " << theta << std::endl;
            Matrix<DType> sub_rot = Matrix<DType>::identity(n);

            Matrix<DType> so2({2,2},
                {cos(theta), -sin(theta),
                sin(theta), cos(theta)});

            sub_rot.setBlock(i-1, i-1, so2);
            // std::cout << mxm::to_string(sub_rot) << std::endl;
            std::cout << "b: \n"<< mxm::to_string(mat)  << std::endl;

            rot = sub_rot.matmul(rot);
            mat = sub_rot.matmul(mat).matmul(sub_rot.T());
        }
    }
    rot = rot.T();
    return ret;
}

} // namespace mxm



#endif // _LINALG_SOLVE_H
