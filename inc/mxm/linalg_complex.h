#if !defined(_LINALG_COMPLEX_H_)
#define _LINALG_COMPLEX_H_

#include <array>
#include <string>
#include <cmath>
#include "common.h"
#include "accessor.h"

namespace mxm
{

template<typename DType, unsigned int N>
class ComplexNumber;

template<typename DType, unsigned int N>
ComplexNumber<DType, N> complexMul(const ComplexNumber<DType, N>& lhs, const ComplexNumber<DType, N>& rhs);

inline std::string complexSymbol(size_t i)
{
    const static std::array<std::string, 4> table{"","i","j","k"};
    return table[i];
}


template<typename DType, unsigned int N>
class ComplexNumber: public LinearData<ComplexNumber<DType, N>, DType>
{
public:
    using ThisType = ComplexNumber<DType,N>;
    ComplexNumber(const std::array<DType, N>& v):data_(v){}
    ComplexNumber(const ThisType& rhs)
    {
        (*this) = rhs;
    }
    ComplexNumber(){}
    ~ComplexNumber(){}

    const DType& operator () (size_t i) const { return data_.at(i); }
    DType& operator () (size_t i) { return const_cast<DType&>(static_cast<const ThisType&>(*this)(i)); }

    constexpr size_t size() const { return N; }

    void operator = (const ThisType& rhs)
    {
        this->traverse([&](size_t i){ (*this)(i) = rhs(i); });
    }

    ThisType& operator *= (DType scalar)  { traverse([&](size_t i){(*this)(i) *= scalar;}); return *this;}
    ThisType operator * (DType scalar) const { return ThisType(*this) *= scalar; }

    ThisType& operator *= (const ThisType& rhs)  { (*this) = complexMul(*this, rhs); return *this;}
    ThisType& operator * (const ThisType& rhs)  { return ThisType(*this) *= rhs;}

    ThisType conjugate() const { ThisType ret(-(*this)); ret(0) *= -1; return ret;}

    std::string str() const
    {
        std::string ret;
        this->traverse([&](size_t i){
            ret += ((*this)(i) >= 0 && i > 0 ? "+" : "");
            ret += (mxm::to_string((*this)(i), 6) + complexSymbol(i)); });
        return ret;
    }

private:
    std::array<DType, N> data_;
};

template<typename DType, unsigned int N>
ComplexNumber<DType, N> complexMul(const ComplexNumber<DType, N>& lhs, const ComplexNumber<DType, N>& rhs)
{
    throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

template<typename DType>
ComplexNumber<DType, 2> complexMul(const ComplexNumber<DType, 2>& lhs, const ComplexNumber<DType, 2>& rhs)
{
    return ComplexNumber<DType, 2>({
        lhs(0) * rhs(0) - lhs(1) * rhs(1),
        lhs(0) * rhs(1) + lhs(1) * rhs(0)});
}

template<typename DType>
ComplexNumber<DType, 4> complexMul(const ComplexNumber<DType, 4>& lhs, const ComplexNumber<DType, 4>& rhs)
{
    const size_t N(4);
    std::array<DType, N> ret_data;
    std::array<DType, N * N> mat = {
        lhs(0), -lhs(1), -lhs(2), -lhs(3),
        lhs(1),  lhs(0), -lhs(3), -lhs(2),
        lhs(2),  lhs(3),  lhs(0), -lhs(1),
        lhs(3), -lhs(2),  lhs(1),  lhs(0)};

    for(size_t i = 0; i < N; i++)
    {
        ret_data[i] = 0;
        for(size_t j = 0; j < N; j++)
        {
            ret_data[i] += mat[i * N + j] * rhs(j);
        }
    }
    return ComplexNumber<DType, 4>(ret_data);
}

using Quaternion = ComplexNumber<FloatType, 4>;
using Complex = ComplexNumber<FloatType, 2>;

} // namespace mxm

#endif // _LINALG_COMPLEX_H_
