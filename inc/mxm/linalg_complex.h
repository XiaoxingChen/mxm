#if !defined(_LINALG_COMPLEX_H_)
#define _LINALG_COMPLEX_H_

#include <array>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include "common.h"
#include "accessor.h"
// #include "linalg_norm.h"
#include <initializer_list>
#include "math_forward_declaration.h"

namespace mxm
{

template<typename DType, unsigned int N>
class Hypercomplex;

template<typename DType, unsigned int N>
Hypercomplex<DType, N> complexMul(const Hypercomplex<DType, N>& lhs, const Hypercomplex<DType, N>& rhs);

inline std::string complexSymbol(size_t i)
{
    const static std::array<std::string, 4> table{"","i","j","k"};
    return table[i];
}

template<typename DType, unsigned int N>
typename Traits<Hypercomplex<DType, N>>::ArithType norm(const Hypercomplex<DType, N>& in);


template<typename DType, unsigned int N>
class Hypercomplex
{
public:
    using ThisType = Hypercomplex<DType,N>;
    Hypercomplex(const std::array<DType, N>& v):data_(v){}
    Hypercomplex(std::initializer_list<DType> v)
    {
        for(size_t i = 0; i < std::distance(v.begin(), v.end()); i++)
        {
            data_.at(i) = *(v.begin() + i);
        }
    }
    Hypercomplex(DType val){ for(auto & d : data_) d = 0; data_.at(0) = val;}
    Hypercomplex(const ThisType& rhs)
    {
        (*this) = rhs;
    }
    Hypercomplex(){}
    ~Hypercomplex(){}

    const DType& operator () (size_t i) const { return data_.at(i); }
    DType& operator () (size_t i) { return const_cast<DType&>(static_cast<const ThisType&>(*this)(i)); }

    static constexpr size_t size() { return N; }

    void traverse(std::function< void(size_t)> f) const { for(size_t i = 0; i < size(); i++) f(i); }
    typename Traits<Hypercomplex<DType, N>>::ArithType norm() const { return mxm::norm(*this); }

    void operator = (const ThisType& rhs)
    {
        this->traverse([&](size_t i){ (*this)(i) = rhs(i); });
    }

    // obj to obj
    ThisType& operator *= (const ThisType& rhs)  { (*this) = complexMul(*this, rhs); return *this;}
    ThisType& operator += (const ThisType& rhs) { this->traverse([&](size_t i) {data_.at(i) += rhs(i);}); return *this; }
    ThisType& operator -= (const ThisType& rhs) { this->traverse([&](size_t i) {data_.at(i) -= rhs(i);}); return *this; }

    ThisType operator * (const ThisType& rhs) const { return ThisType(*this) *= rhs;}
    ThisType operator + (const ThisType& rhs) const { return ThisType(*this) += rhs; }
    ThisType operator - (const ThisType& rhs) const { return ThisType(*this) -= rhs; }


    // obj to scalar, "+/-" only takes effect on real part, "*" take effects on both
    ThisType& operator += (const DType& rhs) { data_.at(0) += rhs; return *this; }
    ThisType& operator -= (const DType& rhs) { data_.at(0) -= rhs; return *this; }

    ThisType operator + (const DType& rhs) const { return ThisType(*this) += rhs; }
    ThisType operator - (const DType& rhs) const { return ThisType(*this) -= rhs; }

    ThisType operator - () const { return ThisType(*this) *= -1; }

    ThisType conj() const { ThisType ret(-(*this)); ret(0) *= -1; return ret;}

    std::string str() const;

private:
    std::array<DType, N> data_;
};

template<typename DType, unsigned int N>
Hypercomplex<DType, N> complexMul(const Hypercomplex<DType, N>& lhs, const Hypercomplex<DType, N>& rhs)
{
    throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template<typename DType>
Hypercomplex<DType, 2> complexMul(const Hypercomplex<DType, 2>& lhs, const Hypercomplex<DType, 2>& rhs)
{
    return Hypercomplex<DType, 2>({
        lhs(0) * rhs(0) - lhs(1) * rhs(1),
        lhs(0) * rhs(1) + lhs(1) * rhs(0)});
}
// Reference:
// [1] Geometric Methods and Applications
// Chapter 9.1
template<typename DType, unsigned int N>
std::enable_if_t<4 == N, std::vector<DType>>
matrixFromLhsComplex(const Hypercomplex<DType, N>& q)
{
    std::vector<DType> ret{
        q(0), -q(1), -q(2), -q(3),
        q(1),  q(0), -q(3), q(2),
        q(2),  q(3),  q(0), -q(1),
        q(3), -q(2),  q(1),  q(0)};
    return ret;
}

// Reference:
// [1] Geometric Methods and Applications
// Chapter 9.1
template<typename DType, unsigned int N>
std::enable_if_t<4 == N, std::vector<DType>>
matrixFromRhsComplex(const Hypercomplex<DType, N>& q)
{
    std::vector<DType> ret{
        q(0),  q(1),  q(2), q(3),
        -q(1), q(0), -q(3), q(2),
        -q(2), q(3),  q(0), -q(1),
        -q(3), -q(2), q(1),  q(0)};
    return ret;
}

template<typename DType>
Hypercomplex<DType, 4> complexMul(const Hypercomplex<DType, 4>& lhs, const Hypercomplex<DType, 4>& rhs)
{
    const size_t N(4);
    std::array<DType, N> ret_data;
    std::vector<DType> mat = matrixFromLhsComplex(lhs);

    for(size_t i = 0; i < N; i++)
    {
        ret_data[i] = 0;
        for(size_t j = 0; j < N; j++)
        {
            ret_data[i] += mat[i * N + j] * rhs(j);
        }
    }
    return Hypercomplex<DType, 4>(ret_data);
}

template<typename DType>
using Quaternion = Hypercomplex<DType, 4>;
template<typename DType>
using Complex = Hypercomplex<DType, 2>;

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value , std::string>::type
to_string(const DType& v, size_t prec);

template<typename DType, unsigned int N>
std::string to_string(const Hypercomplex<DType, N>& v, size_t prec)
{
    std::string ret;
    v.traverse([&](size_t i){
        ret += ((v(i) >= 0 || std::isnan(v(i))) && i > 0 ? "+" : "");
        ret += (mxm::to_string(v(i), prec) + complexSymbol(i)); });
    return ret;
}

template<typename DType, unsigned int N>
std::string to_string(const Hypercomplex<DType, N>& v)
{
    return mxm::to_string(v, 6);
}

template<typename DType, unsigned int N>
std::string Hypercomplex<DType, N>::str() const
{
    return to_string(*this);
}
#if 0
template<typename DType, unsigned int N>
struct NormTraits<Hypercomplex<DType, N>>{using type=DType;};
#endif


template<typename DType, unsigned int N>
Hypercomplex<DType, N> inv(const Hypercomplex<DType, N>& val)
{
    auto norm_v = val.norm();
    return val.conj() * (decltype(norm_v)(1) / (norm_v * norm_v));
}

template<typename DType, unsigned int N>
struct Traits<Hypercomplex<DType, N>>{
    using ArithType = DType;
    using EntryType = DType;
    static Hypercomplex<DType, N> identity() { return Hypercomplex<DType, N>(1); }
    static Hypercomplex<DType, N> zero() { return Hypercomplex<DType, N>(0); }
    // static Hypercomplex<DType, N> inv(const Hypercomplex<DType, N>& val) { return Px<DType>(); }
    // static ArithType norm(const Px<DType>& val) { return ArithType(5); }
    // static std::string to_string(const Px<DType>& val) { return std::string("px"); }
};

template<typename DType, unsigned int N>
typename Traits<Hypercomplex<DType, N>>::ArithType norm(const Hypercomplex<DType, N>& in)
{
    typename Traits<Hypercomplex<DType, N>>::ArithType sum = 0;
    in.traverse([&](auto i) { sum += in(i)*in(i); });
    return sqrt(sum);
}

template<typename DeriveType>
class MatrixBase;


template<typename MatrixType, typename U=void>
struct IsHypercomplexMatrix { static const bool value = false; };
template<typename MatrixType>
struct IsHypercomplexMatrix<MatrixType, std::enable_if_t<Traits<MatrixType>::EntryType::size() >= 0, void>>
{
    static const bool value = std::is_same<
            Hypercomplex<
                typename Traits<MatrixType>::ArithType,
                Traits<MatrixType>::EntryType::size()
            >, typename Traits<MatrixType>::EntryType
        >::value;
};

template<typename MatrixType>
struct IsRealMatrix {
    static const bool value = std::is_arithmetic<typename Traits<MatrixType>::ArithType>::value ;
};

// template<template <class> class MatrixType, class EntryType>
// std::enable_if_t<IsHypercomplexMatrix<MatrixType<EntryType>>::value , Matrix<EntryType>>
// conj(const MatrixBase<MatrixType<EntryType>>& in)
// {
//     const size_t N = Traits<EntryType>::size();
//     Matrix<EntryType> ret(in);
//     ret.traverse([&](auto i, auto j){ret(i,j) = ret(i,j).conj();});
//     return ret;
// }

template<template <class> class MatrixType, class EntryType>
std::enable_if_t<IsHypercomplexMatrix<MatrixType<EntryType>>::value , Matrix<EntryType>>
conj(const MatrixBase<MatrixType<EntryType>>& in)
{
    const size_t N = EntryType::size();
    Matrix<EntryType> ret(in);
    ret.traverse([&](auto i, auto j){ret(i,j) = ret(i,j).conj();});
    return ret;
}

template<typename DType, unsigned int N>
std::enable_if_t<2 == N, DType> imNorm_(const Hypercomplex<DType, N>& in)
{
    return abs(in(1));
}

template<typename DType, unsigned int N>
std::enable_if_t<4 == N, DType> imNorm_(const Hypercomplex<DType, N>& in)
{
    DType im_norm2 = DType(0);
    for(size_t i = 1; i < in.size(); i++)
    {
        im_norm2 += in(i)*in(i);
    }
    return std::sqrt(im_norm2);
}

template<typename DType, unsigned int N>
std::enable_if_t<2 == N, std::array<DType, N-1>> imDir_(const Hypercomplex<DType, N>& in, DType im_norm)
{
    return {1};
}

template<typename DType, unsigned int N>
std::enable_if_t<4 == N, std::array<DType, N-1>> imDir_(const Hypercomplex<DType, N>& in, DType im_norm)
{
    std::array<DType, N-1> ret{};
    if(im_norm < std::numeric_limits<DType>::epsilon())
    {
        ret[0] = 1;
        return ret;
    }
    if(isinf(im_norm))
    {
        for(size_t i = 0; i < ret.size(); i++)
        {
            if(isinf(in(i+1)))
            {
                ret[i] = in(i+1) > 0 ? 1 : -1;
                break;
            }
        }
    }
    for(size_t i = 0; i < ret.size(); i++)
    {
        ret[i] = in(i+1) / im_norm;
    }
    return ret;
}

// template<typename DType, unsigned int N>
// std::enable_if_t<2 == N, std::array<DType, N-1>>
// imDir(const Hypercomplex<DType, N>& in)
// { return {DType(1)}; }

// template<typename DType, unsigned int N>
// std::enable_if_t<4 == N, std::array<DType, N-1>>
// imDir(const Hypercomplex<DType, N>& in)
// {

// }

template<typename DType, unsigned int N>
Hypercomplex<DType, N> log(const Hypercomplex<DType, N>& in)
{
    Hypercomplex<DType, N> ret = Traits<Hypercomplex<DType, N>>::zero();
    DType norm2 = norm(in);
    if(norm2 < std::numeric_limits<DType>::epsilon())
    {
        ret(0) = - std::numeric_limits<DType>::infinity();
        return ret;
    }
    ret(0) = std::log(norm2);
    DType im_norm = imNorm_(in);
    auto im_dir = imDir_(in, im_norm);
    if(im_norm < std::numeric_limits<DType>::epsilon())
    {
        ret(1) = atan2(0, in(0));
        return ret;
    }

    DType angle = std::atan2(im_norm, in(0));

    for(size_t i = 1; i < in.size(); i++)
    {
        ret(i) = in(i) * im_dir[i-1] * angle;
    }
    return ret;
}

template<typename DType, unsigned int N>
Hypercomplex<DType, N> exp(const Hypercomplex<DType, N>& in)
{
    auto ret = Traits<Hypercomplex<DType, N>>::identity();
    ret *= mxm::exp(in(0));
    for(size_t i = 1; i < N; i++)
    {
        auto val = Traits<Hypercomplex<DType, N>>::identity();
        val(0) = mxm::cos(in(i));
        val(i) = mxm::sin(in(i));
        ret *= val;
    }
    return ret;
}

#if 1
template<typename DType, unsigned int N>
Hypercomplex<DType, N> pow(const Hypercomplex<DType, N>& in, DType y)
{
    Hypercomplex<DType, N> ret = Traits<Hypercomplex<DType, N>>::zero();
    DType norm2 = norm(in);
    if(norm2 < std::numeric_limits<DType>::epsilon()) return ret;

    DType im_norm = imNorm_(in);
    if(im_norm < std::numeric_limits<DType>::epsilon())
    {
        ret(0) = std::pow(in(0), y);
        return ret;
    }
    DType norm_out = std::pow(norm(in), y);
    DType angle = std::atan2(imNorm_(in), in(0)) * y;
    auto im_dir = imDir_(in, im_norm);
    ret(0) = std::cos(angle);
    DType k = std::sin(angle) * norm_out;
    for(size_t i = 1; i < in.size(); i++)
        ret(i) = k * im_dir[i-1];
    return ret;
}
#endif

namespace interp
{

// Spherical linear interpolation
// https://en.wikipedia.org/wiki/Slerp
template<typename DType, unsigned int N>
Hypercomplex<DType, N> slerp(const Hypercomplex<DType, N>& q0, const Hypercomplex<DType, N>& q1, DType t)
{
    return pow(q1 * inv(q0), t) * q0;
}

} // namespace interp


} // namespace mxm

#endif // _LINALG_COMPLEX_H_
