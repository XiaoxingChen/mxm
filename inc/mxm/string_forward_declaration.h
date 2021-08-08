#if !defined(_STRING_FORWARD_DECLARATION_H_)
#define _STRING_FORWARD_DECLARATION_H_

#include <type_traits>
#include <string>

namespace mxm
{

// integral to string
template<typename DType>
typename std::enable_if<std::is_integral<DType>::value , std::string>::type
to_string(const DType& v, size_t prec);

// float point to string
template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value , std::string>::type
to_string(const DType& v, size_t prec);

// Hypercomplex to string
template<typename DType, unsigned int N> class Hypercomplex;

template<typename DType, unsigned int N>
std::string to_string(const Hypercomplex<DType, N>& v, size_t prec);

template<typename DType, unsigned int N>
std::string to_string(const Hypercomplex<DType, N>& v);

// Pixel to string
template<typename DType, size_t N>
class PixelType;

template<typename DType, size_t N>
std::string to_string(const PixelType<DType, N>& px, size_t prec);

template<typename DType, size_t N>
std::string to_string(const PixelType<DType, N>& px);

template<class T, size_t N>
std::string to_string(const std::array<T, N>& container, size_t prec=6);

template<uint8_t K>
class BriefDescriptor;

template<uint8_t K>
std::string to_string(const BriefDescriptor<K>& v, size_t prec=6);

} // namespace mxm

#endif // _STRING_FORWARD_DECLARATION_H_
