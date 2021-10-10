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
    return isZero(mat - conj(mat.T()), tol, error) && norm(mat - DType(1)) < tol;
}
} // namespace su

namespace SU
{
template<typename DType>
bool isValid(const Matrix<Complex<DType>>& mat, DType tol=std::numeric_limits<DType>::epsilon(), DType* error=nullptr)
{
    if(! isSquare(mat)) return false;
    return isIdentity(mat.matmul(conj(mat.T())), error, tol);
}

#if 0
template<size_t N, typename DType>
std::enable_if_t<2 == N, Matrix<DType>>
log(const Matrix<Complex<DType>>& mat)
{
    return Matrix<Complex<DType>>({2,2}, {})
}
#endif

} // namespace SU


} // namespace mxm

#endif // _LIE_SPECIAL_UNITARY_H
