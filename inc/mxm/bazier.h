#if !defined(_BAZIER_H_)
#define _BAZIER_H_

#include "linalg.h"

namespace mxm
{
namespace bazier
{
template<typename DType>
inline Matrix<DType> base(size_t n, size_t i_end, const Vector<DType>& t)
{

    Matrix<DType> result({i_end, t.shape(0)});
    for(size_t i = 0; i < i_end; i++)
    {
        size_t coeff = factorial(n) / factorial(i) / factorial(n-i);
        result(Row(i)) = ((-t + 1).pow(i) * t.pow(n-i) * coeff).T();
    }
    return result;
}

template<typename DType>
inline Vector<DType> f(size_t n, const Matrix<DType>& v, const Vector<DType>& t)
{
    if(v.shape(0) != n + 1)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    return v.T().matmul(base(n, v.shape(0), t));
}
} // namespace bazier

} // namespace mxm


#endif // _BAZIER_H_
