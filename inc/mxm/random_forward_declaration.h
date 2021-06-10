#if !defined(_RANDOM_FORWARD_DECLARATION_H_)
#define _RANDOM_FORWARD_DECLARATION_H_

#include <type_traits>
#include <vector>

namespace mxm
{
namespace random
{

// float type
template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, DType>::type
uniform();

// pixel type
template<typename DType, size_t N>
class PixelType;
    
template<typename DType>
typename std::enable_if<
    std::is_same<
        DType, ::mxm::PixelType< typename DType::EntryType , DType::size()>
    >::value, DType
>::type
uniform();

} // namespace random

} // namespace mxm


#endif // _RANDOM_FORWARD_DECLARATION_H_
