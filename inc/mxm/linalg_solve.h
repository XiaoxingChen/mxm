#if !defined(_LINALG_SOLVE_H)
#define _LINALG_SOLVE_H

#include "mxm/linalg_mat_block.h"
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

enum TraverseSeq{
    eUpperTrianglize,
    eUpperHessenbergize,
    eSubdiagonal
};


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

//
// idx_seq: TraverseSeq
// symmetric: output will be QHQ' = mat if enabled.
// output: {Q, R}
template<typename DeriveType>
std::array<Matrix<typename Traits<DeriveType>::EntryType>, 2>
decomposeByRotation(const MatrixBase<DeriveType>& mat_in, TraverseSeq idx_seq=eUpperTrianglize, bool symmetric=false);

template<typename DType>
Matrix<DType> calcMatQFromReflection(const Matrix<DType>& mat);

template<typename DType>
Matrix<DType> solve(const Matrix<DType>& mat_a, const Matrix<DType>& b);

} // namespace qr

// calculate the determinant of a matrix
// available in both real and complex field.
template <typename DeriveType>
typename Traits<DeriveType>::EntryType
det(const MatrixBase<DeriveType>& mat);

// Inversion
template<typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
inv(const MatrixBase<DeriveType>& mat);

// Eigenvalues of 2x2 real matrix
template <template <class> class MatrixType, typename DType>
std::array<Complex<typename Traits<DType>::ArithType>,2> eigvals2x2(const MatrixBase<MatrixType<DType>>& mat);

// A = QTQ', where
// Q is a real, orthogonal matrix,
// T is a real, quasi-upper-triangular matrix that has a block upper-triangular structure
template<typename DType>
Matrix<DType> realSchurDecomposition(const Matrix<DType>& mat, Matrix<DType>* p_orthogonal=nullptr);


// eigenvalue finder for real matrix
template <typename DType>
std::vector<Complex<DType>> eigvals(const Matrix<DType>& mat);

// DType is required to be double or float
// std::is_float_point<DType>::value == true
template <typename DType>
std::array<Matrix<DType>, 2> symmetricEig(const Matrix<DType>& mat);

template<typename DType>
std::array<Matrix<Complex<DType>>, 2>
eig(const Matrix<DType>& mat);

template<typename DeriveType>
std::array<Matrix<typename Traits<DeriveType>::EntryType>, 3>
svd(const MatrixBase<DeriveType>& mat);

} // namespace mxm

#ifdef MXM_HEADER_ONLY
// #include "linalg_solve_inl.h"
#endif

#endif // _LINALG_SOLVE_H
