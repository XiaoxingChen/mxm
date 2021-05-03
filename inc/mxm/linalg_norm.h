#if !defined(_LINALG_NORM_H_)
#define _LINALG_NORM_H_

namespace mxm
{
template<typename T>
struct NormTraits { using type = void; };

// template<T> struct NormTraits<T> { using type = typename std::enable_if<std::is_floating_point<T>::value, T>::type; };
template<> struct NormTraits<float> { using type = float; };
template<> struct NormTraits<double> { using type = double; };

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, DType>::type norm(const DType& in)
{
    return abs(in);
}

#if 0
template<typename InType>
typename NormTraits<InType>::type norm(const InType& in)
{
    throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}
#endif

} // namespace mxm


#endif // _LINALG_NORM_H_
