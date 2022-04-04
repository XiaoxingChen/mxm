// References:
// [1] https://en.wikipedia.org/wiki/Special_unitary_group
// [2] https://www.youtube.com/watch?v=ACZC_XEyg9U Dirac's belt trick, Topology, and Spin 1/2 particles
#if !defined(_LIE_SPECIAL_UNITARY_H)
#define _LIE_SPECIAL_UNITARY_H

#include "linalg.h"

#include <limits>

namespace mxm
{

namespace su
{
template<typename DType>
bool isValid(const Matrix<Complex<DType>>& mat, DType* error=nullptr, DType tol=std::numeric_limits<DType>::epsilon())
{
    if(! isSquare(mat)) return false;
    return isZero(mat + conj(mat.T()), tol, error);
}
#if 1
// Reference
// https://en.wikipedia.org/wiki/Special_unitary_group
template<size_t N, typename DeriveType>
std::enable_if_t<2 == N, Matrix<typename Traits<DeriveType>::EntryType>>
exp(const MatrixBase<DeriveType>& mat)
{
    using DType = typename Traits<DeriveType>::ArithType;
    DType x = ((mat(0,1) + mat(1,0)) * 0.5)(1);
    DType y = ((mat(1,0) - mat(0,1)) * 0.5)(0);
    DType z = ((mat(0,0) - mat(1,1)) * 0.5)(1);
    DType half_angle = sqrt(x*x  + y*y + z*z);
    if(half_angle < eps<DType>())
        return Matrix<typename Traits<DeriveType>::EntryType>::identity(N);

    DType inv_half_angle = DType(1) / half_angle;
    DType s = sin(half_angle);
    DType c = cos(half_angle);

    Matrix<Complex<DType>> uw({2,2}, {1, 0, 0, 1}, ROW);
    Matrix<Complex<DType>> ux({2,2}, {{0,1}, 0, 0, {0,-1}}, ROW);
    Matrix<Complex<DType>> uy({2,2}, {0, 1, -1, 0}, ROW);
    Matrix<Complex<DType>> uz({2,2}, {0, {0,1}, {0,1}, 0}, ROW);

    return uw * c + (ux * x + uy * y + uz * z) * s * inv_half_angle;
}
#endif

} // namespace su

namespace SU
{
template<typename DType>
bool isValid(const Matrix<Complex<DType>>& mat, DType* error=nullptr, DType tol=std::numeric_limits<DType>::epsilon())
{
    if(! isSquare(mat)) return false;
    return isIdentity(mat.matmul(conj(mat.T())), error, tol) && norm(mxm::det(mat) - DType(1)) < tol;
}

#if 1
// Reference
// https://en.wikipedia.org/wiki/Special_unitary_group
template<size_t N, typename DeriveType>
std::enable_if_t<2 == N, Matrix<typename Traits<DeriveType>::EntryType>>
log(const MatrixBase<DeriveType>& mat)
{
    using DType = typename Traits<DeriveType>::ArithType;
    DType w = ((mat(0,0) + mat(1,1)) * 0.5)(0);
    DType x = ((mat(0,0) - mat(1,1)) * 0.5)(1);
    DType y = ((mat(0,1) - mat(1,0)) * 0.5)(0);
    DType z = ((mat(0,1) + mat(1,0)) * 0.5)(1);
    DType norm_imagine = sqrt(x*x + y*y + z*z);
    if(norm_imagine < eps<DType>()) return Matrix<Complex<DType>>::zeros({2,2});
    DType half_angle = atan2(norm_imagine, w);
    DType inv_sin_half_angle = DType(1.) / norm_imagine;
    Matrix<Complex<DType>> u1({2,2}, {0, {0,1}, {0,1}, 0}, ROW);
    Matrix<Complex<DType>> u2({2,2}, {0, -1, 1, 0}, ROW);
    Matrix<Complex<DType>> u3({2,2}, {{0,1}, 0, 0, {0,-1}}, ROW);
    return (u1 * x + u2 * y + u3 * z) * inv_sin_half_angle * half_angle;
}
#endif


template<typename DType>
Matrix<DType> inv(const Matrix<DType>& mat)
{
    return conj(mat.T());
}

} // namespace SU


} // namespace mxm

#endif // _LIE_SPECIAL_UNITARY_H
