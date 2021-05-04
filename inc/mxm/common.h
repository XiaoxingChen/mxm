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

namespace mxm{
using FloatType = float;

inline constexpr FloatType eps() {return 1e-7;}

inline constexpr FloatType tMin() {return 1e-4;}
inline constexpr FloatType tMax() {return 10000.;}

// template<typename DType>

template<typename DType>
inline typename std::enable_if<std::is_integral<DType>::value , std::string>::type
to_string(const DType& v, size_t prec=6)
{
    return std::to_string(v);
}

template<typename DType>
inline typename std::enable_if<std::is_floating_point<DType>::value , std::string>::type
to_string(const DType& v, size_t prec=6)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(prec) << v;
    return stream.str();
}

template<template <class, class> class Container, class T, class Alloc>
inline std::string to_string(const Container<T, Alloc>& container, size_t prec=6)
{
    std::string ret;
    for(size_t i = 0; i < container.size(); i++)
    {
        ret += to_string(container.at(i), prec);
        ret += (i == container.size() - 1 ? "" : ret.back() == '\n' ? "\n" : " ");
    }
    return ret;
}

template<class T, size_t N>
inline std::string to_string(const std::array<T, N>& container, size_t prec=6)
{
    std::string ret;
    for(size_t i = 0; i < container.size(); i++)
    {
        ret += to_string(container.at(i), prec);
        ret += (i == container.size() - 1 ? "" : ret.back() == '\n' ? "\n" : " ");
    }
    return ret;
}



}//namespace mxm


#endif // _COMMON_H_
