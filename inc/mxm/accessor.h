#if !defined(_ACCESSOR_H_)
#define _ACCESSOR_H_

#include <functional>
#include <cmath>
#include <string>
#include "common.h"

namespace mxm
{

using LoopN1 = std::function< void(size_t)>;
using LoopN2 = std::function< void(size_t, size_t)>;


// template<typename DType>
// inline std::string to_string(const DType& v, size_t prec=6)
// {
//     std::stringstream stream;
//     stream << std::fixed << std::setprecision(prec) << v;
//     return stream.str();
// }

//DType can only be float, double or int
template<typename DeriveType, typename DType>
class LinearData
{
public:
    void traverse(std::function< void(size_t)> f) const { for(size_t i = 0; i < deriveThis()->size(); i++) f(i); }

    DeriveType* deriveThis() {return reinterpret_cast<DeriveType*>(this);}
    const DeriveType* deriveThis() const {return reinterpret_cast<const DeriveType*>(this);}
    const DType& at(size_t i) const {return (*deriveThis())(i);}
    DType& at(size_t i) {return (*deriveThis())(i);}

    DType norm() const { DType sum; traverse([&, this](size_t i) { DType v = at(i); sum += v * v;}); return sqrt(sum);}
    DeriveType& normalize() { return (*deriveThis()) *= (1./DeriveType::norm()); }
    DeriveType normalized() const { return DeriveType(*deriveThis()).normalize(); }
    std::string str() const { std::string ret; traverse([&](size_t i) {ret += to_string(at(i)) + " ";}); return ret; }

    // obj to obj
    DeriveType& operator += (const DeriveType& rhs) { traverse([&](size_t i) {at(i) += rhs(i);}); return *deriveThis(); }
    DeriveType& operator -= (const DeriveType& rhs) { traverse([&](size_t i) {at(i) -= rhs(i);}); return *deriveThis(); }
    DeriveType& operator *= (const DeriveType& rhs) { traverse([&](size_t i) {at(i) *= rhs(i);}); return *deriveThis(); }

    DeriveType operator + (const DeriveType& rhs) const { return DeriveType(*deriveThis()) += rhs; }
    DeriveType operator - (const DeriveType& rhs) const { return DeriveType(*deriveThis()) -= rhs; }
    DeriveType operator * (const DeriveType& rhs) const { return DeriveType(*deriveThis()) *= rhs; }

    DeriveType operator - () { return DeriveType(*deriveThis()) *= -1; }

    // obj to scalar
    DeriveType& operator += (const DType& rhs) { traverse([&](size_t i) {at(i) += rhs;}); return *deriveThis(); }
    DeriveType& operator -= (const DType& rhs) { traverse([&](size_t i) {at(i) -= rhs;}); return *deriveThis(); }
    DeriveType& operator *= (const DType& rhs) { traverse([&](size_t i) {at(i) *= rhs;}); return *deriveThis(); }

    DeriveType operator + (const DType& rhs) const { return DeriveType(*deriveThis()) += rhs; }
    DeriveType operator - (const DType& rhs) const { return DeriveType(*deriveThis()) -= rhs; }
    DeriveType operator * (const DType& rhs) const { return DeriveType(*deriveThis()) *= rhs; }

    // (13.5.3 Assignment)
    // An assignment operator shall be implemented by a non-static member function with exactly one parameter.
    // Because a copy assignment operator operator= is implicitly declared for a a class if not declared by the user,
    // a base class assignment operator is always hidden by the copy assignment operator of the derived class.
    void operator = (const DeriveType& rhs) { /*none sense*/ }

    bool operator == (const DeriveType& rhs) const
    {
        bool ret(true);
        traverse([&](size_t i) { ret &= (at(i) == rhs(i));});
        return ret;
    }

    bool operator != (const DeriveType& rhs) const {return !((*this) == rhs);}


};

template<typename DeriveType, typename DType> inline DeriveType operator + (const DType& scalar, const LinearData<DeriveType, DType>& rhs ) {return rhs + scalar;}
template<typename DeriveType, typename DType> inline DeriveType operator * (const DType& scalar, const LinearData<DeriveType, DType>& rhs ) {return rhs * scalar;}

} // namespace mxm


#endif // _ACCESSOR_H_
