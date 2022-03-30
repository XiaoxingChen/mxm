#if !defined(_COMMON_H_)
#define _COMMON_H_

#ifdef MXM_COMPILED_LIB
#undef MXM_HEADER_ONLY
#if defined(_WIN32) && defined(MXM_SHARED_LIB)
#ifdef MXM_EXPORTS
#define MXM_API __declspec(dllexport)
#else
#define MXM_API __declspec(dllimport)
#endif
#else // !defined(_WIN32) || !defined(MXM_SHARED_LIB)
#define MXM_API
#endif
#define MXM_INLINE
#else // !defined(MXM_COMPILED_LIB)
#define MXM_API
#define MXM_HEADER_ONLY
#define MXM_INLINE inline
#endif // #ifdef MXM_COMPILED_LIB

#include <limits>

namespace mxm{
using FloatType = float;

template<typename DType=float>
inline constexpr DType eps() {return std::numeric_limits<DType>::epsilon();}

inline constexpr FloatType tMin() {return 1e-4;}
inline constexpr FloatType tMax() {return 10000.;}

// template<typename DType>


}//namespace mxm


#endif // _COMMON_H_
