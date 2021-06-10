#if !defined(_LINALG_NORM_H_)
#define _LINALG_NORM_H_
#include <type_traits>
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

template<typename T, typename>
struct Traits;

template<typename T>
struct MatrixBase;

template<typename DType,
    typename=std::enable_if_t<std::is_floating_point<DType>::value, void>>
typename Traits<DType>::ArithType norm(const DType& val){ return std::abs(val); }

template<typename DeriveType, typename>
typename Traits<DeriveType>::ArithType
norm(const MatrixBase<DeriveType>& mat)
{
    using EntryType = typename Traits<DeriveType>::EntryType;
    using ArithType = typename Traits<DeriveType>::ArithType;
    auto& self = reinterpret_cast<const DeriveType&>(mat);
    typename Traits<DeriveType>::ArithType sum2(0);
    self.traverse([&](size_t i, size_t j){
        auto n = norm(self(i,j)); sum2 += n*n;
    });
    return sqrt(sum2);
}

} // namespace mxm


#endif // _LINALG_NORM_H_
