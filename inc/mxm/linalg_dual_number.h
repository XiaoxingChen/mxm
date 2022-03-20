#if !defined(__LINALG_DUAL_NUMBER_H__)
#define __LINALG_DUAL_NUMBER_H__

// References:
// [1] https://www.youtube.com/watch?v=ceaNqdHdqtg
// [2] https://en.wikipedia.org/wiki/Dual_number

#include <array>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include "common.h"
// #include "linalg_norm.h"
#include <initializer_list>

namespace mxm
{
template <typename DType> class DualNumber;

template<typename DType>
DualNumber<DType> dualNumberMul(const DualNumber<DType>& lhs, const DualNumber<DType>& rhs)
{
    return DualNumber<DType>{lhs(0) * rhs(0), lhs(0) * rhs(1) + lhs(1) * rhs(0)};
}

template<typename DType>
DualNumber<DType> inv(const DualNumber<DType>& val)
{
    DType denominator = val.conj() * val;
    return val.conj() * (DType(1.) / denominator);
}

template<typename DType>
struct Traits<DualNumber<DType>>{
    using ArithType = DType;
    using EntryType = DType;
    static DualNumber<DType> identity() { return DualNumber<DType>(1); }
    static DualNumber<DType> zero() { return DualNumber<DType>(0); }
};

#if 0
template<typename DType>
bool isZero(
    const DualNumber<DType>& val,
    typename Traits<DType>::ArithType* p_error,
    typename Traits<DType>::ArithType tol)
{
    using ArithType = typename Traits<DType>::ArithType;
    ArithType error = std::sqrt(val(0)*val(0) + val(1) * val(1));
    if(p_error) *p_error = error;
    return error < tol;
}
#endif

template<typename DType>
std::string to_string(const DualNumber<DType>& v, size_t prec)
{
    std::string ret;
    v.traverse([&](size_t i){
        ret += ((v(i) >= 0 || std::isnan(v(i))) && i > 0 ? "+" : "");
        ret += (mxm::to_string(v(i), prec) + "eps"); });
    return ret;
}

template <typename DType>
class DualNumber
{
public:
    using ThisType = DualNumber<DType>;
    DualNumber(std::initializer_list<DType> v)
    {
        size_t list_len = std::distance(v.begin(), v.end());
        for(size_t i = 0; i < data_.size() ; i++)
        {
            data_.at(i) = i < list_len ? *(v.begin() + i) : DType(0);
        }
    }
    DualNumber(DType val){ for(auto & d : data_) d = 0; data_.at(0) = val;}
    DualNumber(const ThisType& rhs)
    {
        (*this) = rhs;
    }
    DualNumber(){}
    ~DualNumber(){}

    const DType& operator () (size_t i) const { return data_.at(i); }
    DType& operator () (size_t i) { return const_cast<DType&>(static_cast<const ThisType&>(*this)(i)); }

    // obj to obj
    ThisType& operator *= (const ThisType& rhs)  { (*this) = dualNumberMul(*this, rhs); return *this;}
    ThisType& operator /= (const ThisType& rhs)  { (*this) *= inv(rhs); return *this;}
    ThisType& operator += (const ThisType& rhs) { this->traverse([&](size_t i) {data_.at(i) += rhs(i);}); return *this; }
    ThisType& operator -= (const ThisType& rhs) { this->traverse([&](size_t i) {data_.at(i) -= rhs(i);}); return *this; }

    ThisType operator * (const ThisType& rhs) const { return ThisType(*this) *= rhs;}
    ThisType operator / (const ThisType& rhs) const { return ThisType(*this) /= rhs;}
    ThisType operator + (const ThisType& rhs) const { return ThisType(*this) += rhs; }
    ThisType operator - (const ThisType& rhs) const { return ThisType(*this) -= rhs; }


    // obj to scalar, "+/-" only takes effect on real part, "*" take effects on both
    ThisType& operator += (const DType& rhs) { data_.at(0) += rhs; return *this; }
    ThisType& operator -= (const DType& rhs) { data_.at(0) -= rhs; return *this; }
    ThisType& operator *= (const DType& rhs) { this->traverse([&](size_t i) {data_.at(i) *= rhs;}); return *this; }
    ThisType& operator /= (const DType& rhs) { this->traverse([&](size_t i) {data_.at(i) /= rhs;}); return *this; }

    ThisType operator + (const DType& rhs) const { return ThisType(*this) += rhs; }
    ThisType operator - (const DType& rhs) const { return ThisType(*this) -= rhs; }
    ThisType operator * (const DType& rhs) const { return ThisType(*this) *= rhs; }
    ThisType operator / (const DType& rhs) const { return ThisType(*this) /= rhs; }

    ThisType operator - () const { return ThisType(*this) *= -1; }

    ThisType conj() const { return ThisType{data_.at(0), data_.at(1) * DType(-1)};}

    std::string str() const;

    void traverse(std::function< void(size_t)> f) const { for(size_t i = 0; i < size(); i++) f(i); }


    template <typename T> operator
    T () const { return static_cast<T>(data_.at(0)); }

private:
    static constexpr size_t size() { return 2; }


    std::array<DType, 2> data_;
};

template<typename DType> DualNumber<DType> operator + (DType lhs, const DualNumber<DType>& rhs) { return rhs + lhs; }
template<typename DType> DualNumber<DType> operator - (DType lhs, const DualNumber<DType>& rhs) { return (-rhs) + lhs; }
template<typename DType> DualNumber<DType> operator * (DType lhs, const DualNumber<DType>& rhs) { return rhs * lhs; }
template<typename DType> DualNumber<DType> operator / (DType lhs, const DualNumber<DType>& rhs) { return inv(rhs) * lhs; }

template<typename DType> DualNumber<DType>
sin(const DualNumber<DType>& val)
{
    return {std::sin(val(0)), std::cos(val(0)) * val(1)};
}

template<typename DType> DualNumber<DType>
cos(const DualNumber<DType>& val)
{
    return DualNumber<DType>{std::cos(val(0)), std::sin(val(0)) * val(1) * DType(-1)};
}

template<typename DType> DType
abs(const DualNumber<DType>& val)
{
    return std::abs(val(0));
}

template<typename DType> DualNumber<DType>
sqrt(const DualNumber<DType>& val)
{
    DType val_0_sqrt = std::sqrt(val(0));
    if(val_0_sqrt < eps<DType>())
        return {0, DType(0.5) * val(1) * INFINITY};
    return {val_0_sqrt, DType(0.5) * val(1) / val_0_sqrt};
}

template<typename DType> DualNumber<DType>
asin(const DualNumber<DType>& val)
{
    return std::asin(val(0)) + DType(1.) / mxm::sqrt(1 - val * val) * DualNumber<DType>({0, val(1)});
}

template<typename DType> DualNumber<DType>
acos(const DualNumber<DType>& val)
{
    return std::acos(val(0)) + DType(-1.) / mxm::sqrt(1 - val * val) * DualNumber<DType>({0, val(1)});
}

template<typename DType> DualNumber<DType>
pow(const DualNumber<DType>& val, DType exp)
{
    return {
        std::pow<DType>(val(0), exp),
        exp * std::pow<DType>(val(0), exp - DType(1)) * val(1)};
}

template<typename DType> DualNumber<DType>
exp(const DualNumber<DType>& val)
{
    DType exp_a = std::exp<DType>(val(0));
    return {exp_a,  exp_a * val(1)};
}

template<typename DType> DualNumber<DType>
log(const DualNumber<DType>& val)
{
    return {val(0),  val(1) / val(0)};
}



} // namespace mxm


#endif // __LINALG_DUAL_NUMBER_H__
