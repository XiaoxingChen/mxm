#if !defined(_LINALG_NORM_H_)
#define _LINALG_NORM_H_

namespace mxm
{
template<typename T>
struct NormTraits { using type = void; };

template<> struct NormTraits<double> { using type = double; };
template<> struct NormTraits<float> { using type = float; };

#if 0
template<typename InType>
typename NormTraits<InType>::type norm(const InType& in)
{
    throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}
#endif

} // namespace mxm


#endif // _LINALG_NORM_H_
