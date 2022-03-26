#if !defined(_MXM_STRING_H_)
#define _MXM_STRING_H_

#include <string>
#include <array>
#include "string_forward_declaration.h"

// #include "cv_pixel.h"

namespace mxm
{

template<typename DType>
inline typename std::enable_if<std::is_integral<DType>::value , std::string>::type
to_string(const DType& v, size_t prec)
{
    return std::to_string(v);
}

template<typename DType>
inline typename std::enable_if<std::is_floating_point<DType>::value , std::string>::type
to_string(const DType& v, size_t prec)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(prec) << v;
    return stream.str();
}

template<template <class, class> class Container, class T, class Alloc>
std::string to_string(const Container<T, Alloc>& container, size_t prec=6)
{
    std::string ret;
    for(auto it = container.begin(); it != container.end(); it++)
    {
        ret += to_string(*it, prec);
        ret += (it == container.end() - 1 ? "" : ret.back() == '\n' ? "\n" : " ");
    }
    return ret;
}

template<template <class, class, class> class Container, class T, class Comp, class Alloc>
std::string to_string(const Container<T, Comp, Alloc>& container, size_t prec)
{
    std::string ret;
    for(auto it = container.begin(); it != container.end(); it++)
    {
        ret += to_string(*it, prec);
        ret += (it == std::prev(container.end()) ? "" : ret.back() == '\n' ? "\n" : " ");
    }
    return ret;
}

template<class T, size_t N>
std::string to_string(const std::array<T, N>& container, size_t prec)
{
    std::string ret;
    for(size_t i = 0; i < container.size(); i++)
    {
        ret += to_string(container.at(i), prec);
        ret += (i == container.size() - 1 ? "" : ret.back() == '\n' ? "\n" : " ");
    }
    return ret;
}
#if 0
template<typename DType>
std::string to_string(const Matrix<DType>& mat, size_t prec)
{
    std::string ret;
    mat.traverse([&](size_t i, size_t j){
        ret += (mxm::to_string(mat(i,j), prec) + (j == mat.shape(1) - 1 ? "\n" : " ")); });
    return ret;
}
#elif 1
template<typename DeriveType>
std::string to_string(const MatrixBase<DeriveType>& mat_in, size_t prec=6)
{
    auto & mat = reinterpret_cast<const DeriveType&>(mat_in);
    std::string ret;
    mat.traverse([&](size_t i, size_t j){
        ret += (mxm::to_string(mat(i,j), prec) + (j == mat.shape(1) - 1 ? "\n" : " "));
        // ret += (mxm::to_string(mat(i,j), prec) + (j == mat.shape(1) - 1 ? ",\n" : ", "));
    });
    return ret;
}
#endif
} // namespace mxm


#endif // _MXM_STRING_H_
