#if !defined(_LINALG_VEC_H)
#define _LINALG_VEC_H

#include <vector>
#include <array>
#include <functional>
#include <numeric>
#include <array>
#include <string>
#include <assert.h>
#include <math.h>
#include "mxm/common.h"
#include "mxm/linalg_mat.h"


namespace mxm
{

// class Pixel;

template<typename DType>
class Vector: public Matrix<DType>
{
public:
    using ThisType = Vector<DType>;
    using BaseType = Matrix<DType>;

    Vector(): BaseType(){}
    Vector(size_t size): BaseType({size, 1}){}
    Vector(const std::vector<DType>& v): BaseType({v.size(), 1}, v){}

    Vector(const BaseType& mat)
        : BaseType(mat.shape(1) == 1 ? mat : mat.T())
    {
        if(mat.shape(0) != 1 && mat.shape(1) != 1)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    size_t size() const { return BaseType::shape(0); }
    DType sum() const
    {
        DType sum(0);
        for(int i = 0; i < size(); i++) sum += (*this)(i);
        return sum;
    }


    virtual DType& operator () (size_t i) { return BaseType::operator()(i, 0); }
    virtual const DType& operator () (size_t i) const { return BaseType::operator()(i, 0); }

    DType dot(const ThisType& rhs) const
    {
        if(size() != rhs.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        DType sum(0);
        for(int i = 0; i < size(); i++) sum += (*this)(i) * rhs(i);
        return sum;
    }

    static ThisType zeros(size_t n) { return ThisType(n); }
    static ThisType ones(size_t n) { return ThisType(n) + 1; }

    operator std::vector<DType> () const
    {
        std::vector<DType> ret;
        for(size_t i = 0; i < size(); i++) ret.push_back((*this)(i));
        return ret;
    }
};

using Vec = Vector<FloatType>;
#if 0
class UnitVec: public Vec
{
public:
    UnitVec(size_t size): Vec(size) { Vec::operator()(0) = 1;}
    UnitVec(const Vec& v): Vec(v) {this->normalize();}
    UnitVec(const Mat& m): Vec(m) {this->normalize();}
    UnitVec(const std::vector<FloatType>& v): Vec(v) {this->normalize();}

    Mat& operator *= (FloatType scalar) = delete;
    Mat& operator += (FloatType scalar) = delete;
    Mat& operator -= (FloatType scalar) = delete;
    Mat& operator *= (const Mat& rhs) = delete;
    Mat& operator += (const Mat& rhs) = delete;
    Mat& operator -= (const Mat& rhs) = delete;

    const FloatType& operator () (size_t i) const { return Vec::operator()(i); }

    // FloatType& operator () (size_t i) {}
    // FloatType& operator () (size_t i, size_t j) {}

private:
    using Mat::normalized;
    using Mat::normalize;
};

using VecIn = const Vec&;
using UnitVecIn = const UnitVec&;
#endif

} // namespace mxm


#endif // _LINALG_VEC_H
