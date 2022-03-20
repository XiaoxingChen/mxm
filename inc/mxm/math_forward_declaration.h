#if !defined(__MATH_FORWARD_DECLARATION_H__)
#define __MATH_FORWARD_DECLARATION_H__

#include <type_traits>

namespace mxm
{

template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> sin(DType val);
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> cos(DType val);
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> abs(DType val);
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> sqrt(DType val);
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> pow(DType val, DType exp);
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> exp(DType val);

template<typename DType>
std::enable_if_t<std::is_floating_point<DType>::value, bool>
isZero(DType val, DType* p_error, DType tol);

template <typename DType> class DualNumber;

template<typename DType> DualNumber<DType> sin(const DualNumber<DType>& val);
template<typename DType> DualNumber<DType> cos(const DualNumber<DType>& val);
template<typename DType> DType abs(const DualNumber<DType>& val);
template<typename DType> DualNumber<DType> sqrt(const DualNumber<DType>& val);
template<typename DType> DualNumber<DType> pow(const DualNumber<DType>& val, DType exp);
template<typename DType> DualNumber<DType> exp(const DualNumber<DType>& val);
template<typename DType> DualNumber<DType> log(const DualNumber<DType>& val);

#if 0
template<typename DType>
bool isZero(
    const DualNumber<DType>& val,
    typename Traits<DType>::ArithType* p_error=nullptr,
    typename Traits<DType>::ArithType tol=eps<typename Traits<DType>::ArithType>());
#endif
} // namespace mxm


#endif // __MATH_FORWARD_DECLARATION_H__
