#if !defined(_LINALG_NORM_H_)
#define _LINALG_NORM_H_

namespace mxm
{
#if 0
template<typename T, typename=void>
struct NormTraits { using type = void; };

template<typename T> struct NormTraits<T, std::enable_if_t<std::is_arithmetic<T>::value, void>> { using type = T; };

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, DType>::type norm(const DType& in)
{
    return abs(in);
}


template<typename InType>
typename NormTraits<InType>::type norm(const InType& in)
{
    throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}


template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, DType>::type
inv(const DType& val) { return DType(1) / val; }
#endif

} // namespace mxm


#endif // _LINALG_NORM_H_
