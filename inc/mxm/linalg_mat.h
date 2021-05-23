#if !defined(_LINALG_MAT_H)
#define _LINALG_MAT_H

#include <vector>
#include <functional>
#include <numeric>
#include <array>
#include <string>
#include <assert.h>
#include <math.h>
#include <sstream>
#include <iomanip>

#include "mxm/common.h"
#include "linalg_complex.h"



namespace mxm
{

size_t indexConvert2D(size_t i, size_t j, bool major, size_t shape_i, size_t shape_j);

template<typename DType> class Matrix;

template<typename DType>
std::string to_string(const Matrix<DType>& mat, size_t prec=6);

class AutoShape
{
public:
    AutoShape(const std::array<size_t, 2>& shape) { for(auto & i : {0,1}) raw_shape_[i] = static_cast<int64_t>(shape[i]); }
    AutoShape(std::initializer_list<size_t> shape) { for(auto & i : {0,1}) raw_shape_[i] = static_cast<int64_t>(*(shape.begin() + i)); }
    AutoShape(int64_t row, int64_t col): raw_shape_{row, col} {}
    std::array<size_t, 2> deduct(size_t total_num) const;

private:
    std::array<int64_t, 2> raw_shape_;
};
AutoShape fixRow(size_t n);
AutoShape fixCol(size_t n);

class Block;
template<typename DType> class MatrixRef;
// using MatRef = MatrixRef<FloatType>;

using Shape = std::array<size_t, 2>;
template<typename DType>
class Matrix
{
public:

    static const bool ROW = 0;
    static const bool COL = 1;
    using ThisType = Matrix<DType>;
    using DataPtr = std::vector<DType>*;
    using EntryType = DType;

    //
    // interfaces
    size_t shape(uint8_t i) const {return shape_.at(i);}
    const Shape& shape() const {return shape_;}
    bool square() const {return shape(0) == shape(1);}
    bool majorAxis() const {return major_;}
    bool& majorAxis() {return major_;}
    virtual ThisType& owner() { return *this; }
    virtual const ThisType& owner() const { return *this; }
    virtual const Shape& absOffset() const { static Shape zero({0,0});  return zero; }

    const DType* data() { return data_.data(); }

    //
    // constructors
    Matrix():shape_({0,0}), data_(), major_(ROW) {}

    Matrix(const AutoShape& _shape, const std::vector<DType>& data={}, bool major=ROW)
        :shape_(_shape.deduct(data.size())), data_(shape_[0] * shape_[1]), major_(major)
    {
        if(data.size() > 0)
            data_ = data;
    }

    Matrix(const AutoShape& _shape, std::vector<DType>&& data, bool major=ROW)
        :shape_(_shape.deduct(data.size())), data_(std::move(data)), major_(major) {}

    Matrix(const ThisType& rhs)
        :shape_(rhs.shape()), data_(rhs.shape(0) * rhs.shape(1)), major_(rhs.majorAxis())
    {
        (*this) = rhs;
    }

    template<typename T>
    Matrix(const Matrix<T>& rhs)
        :shape_(rhs.shape()), data_(rhs.shape(0) * rhs.shape(1)), major_(rhs.majorAxis())
    {
        traverse([&](auto i, auto j){
            (*this)(i,j) = DType(rhs(i,j));
        });
    }

    Matrix(ThisType&& rhs)
        :shape_({0,0}), data_(0), major_(ROW)
    {
        (*this) = std::move(rhs);
    }

    virtual void operator = (ThisType&& rhs)
    {
        if(&rhs.owner() == &rhs) // type(rhs) != MatrixRef<DType>
        {
            shape_= rhs.shape_;
            data_.swap(rhs.data_);
            major_ = rhs.major_;

            rhs.shape_ = {0,0};
            return;
        }

    #if 0
        // disable move construct from MatrixRef. See testMatRef()-test09
        if(rhs.owner().shape(0) * rhs.owner().shape(1) == rhs.shape(0) * rhs.shape(1)) // rhs full shape MatrixRef<DType>
        {
            shape_= rhs.shape_;
            data_.swap(rhs.owner().data_);
            major_ = rhs.major_;

            rhs.shape_ = {0,0};
            rhs.owner().shape_ = {0,0};
            return ;
        }
    #endif

        (*this) = rhs;
    }

    virtual void operator = (const ThisType& rhs)
    {
        if(&rhs == this) return; // deal with `mat = mat;`

        if(&rhs.owner() == this) // self assignment
        {
            if(rhs.majorAxis() != majorAxis()
                && rhs.shape() == Shape({shape_[1], shape_[0]}))
            {   //deal with `mat = mat.T();`
                major_ = ! major_;
            }else if(rhs.majorAxis() == majorAxis() && rhs.shape() == shape())
            {   // deal with `mat = mat(Block({0,end()},{0,end()}));
                return; //pass
            }else
            {   // deal with `mat = mat(Block({1,end()}, {1,end()})).T();`
                ThisType tmp(rhs);
                (*this) = std::move(tmp);
            }
        }else // rhs != this
        {
            shape_ = rhs.shape();
            data_.resize(shape(0) * shape(1));
            major_ = rhs.majorAxis();
            traverse([&](size_t i, size_t j){ (*this)(i,j) = rhs(i,j); });
        }
    }

    //
    // static methods
    //
    static ThisType zeros(const Shape& _shape) { return ThisType(_shape) *= DType(0); }
    static ThisType ones(const Shape& _shape) { return zeros(_shape) += DType(1); }
    static ThisType identity(size_t n)
    {
        ThisType mat = zeros({n,n});
        for(int i = 0; i < n; i++) mat(i,i) += DType(1);
        return mat;
    }

    //
    // basic accessor
    //
    virtual const DType& operator () (size_t i, size_t j) const
    {
        return data_.at(
            indexConvert2D(i,j,major_, shape_[0], shape_[1]));
    }

    virtual DType& operator () (size_t i, size_t j)
    {
        return data_.at(
            indexConvert2D(i,j,major_, shape_[0], shape_[1]));
    }

    template<typename Op>
    void traverse(Op f) const
    {
        for(size_t i = 0; i < shape_[0]; i++)
            for(size_t j = 0; j < shape_[1]; j++)
                f(i, j);
    }

    //
    // arithmetic operators
    //
    ThisType& operator *= (DType scalar)  { traverse([&](size_t i, size_t j){(*this)(i,j) *= scalar;}); return *this;}
    ThisType& operator += (DType scalar)  { traverse([&](size_t i, size_t j){(*this)(i,j) += scalar;}); return *this;}
    ThisType& operator -= (DType scalar)  { traverse([&](size_t i, size_t j){(*this)(i,j) -= scalar;}); return *this;}

    ThisType operator * (DType scalar) const { return ThisType(*this) *= scalar; }
    ThisType operator + (DType scalar) const { return ThisType(*this) += scalar; }
    ThisType operator - (DType scalar) const { return ThisType(*this) -= scalar; }

    ThisType& operator *= (const ThisType& rhs)  { traverse([&](size_t i, size_t j){(*this)(i,j) *= rhs(i,j);}); return *this;}
    ThisType& operator += (const ThisType& rhs)  { traverse([&](size_t i, size_t j){(*this)(i,j) += rhs(i,j);}); return *this;}
    ThisType& operator -= (const ThisType& rhs)  { traverse([&](size_t i, size_t j){(*this)(i,j) -= rhs(i,j);}); return *this;}

    ThisType operator * (const ThisType& rhs) const { return ThisType(*this) *= rhs; }
    ThisType operator + (const ThisType& rhs) const { return ThisType(*this) += rhs; }
    ThisType operator - (const ThisType& rhs) const { return ThisType(*this) -= rhs; }

    ThisType operator -() const        { return ThisType(*this) *= -1;}

    bool operator == (const ThisType& rhs) { bool ret(true); traverse([&](size_t i, size_t j){ret = (ret && (*this)(i,j) == rhs(i,j));}); return ret;}
    bool operator != (const ThisType& rhs) { return ! (*this == rhs); }

#if 1
    template<typename DTypeRHS>
    Matrix<decltype(DType() * DTypeRHS())> operator * (const Matrix<DTypeRHS>& rhs) const
    {
        Matrix<decltype(DType() * DTypeRHS())> ret(shape());
        ret.traverse([&](size_t i, size_t j){ret(i,j) = (*this)(i,j) * rhs(i,j);});
        return ret;
    }
#endif

    DType det() const;
    ThisType inv() const;
    // Matrix<Complex<DType>> eigvals() const;

    DType trace() const
    {
        if(!square()) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        DType ret(0);
        for(size_t i = 0; i < shape(0); i++) ret += (*this)(i,i);
        return ret;
    }

    typename NormTraits<DType>::type norm(uint8_t p = 2) const;

    const ThisType& normalize(uint8_t p = 2)
    {
        auto n = norm(p);
        if(n < eps()*eps()) { return *this; }
        return ((*this) *= mxm::inv(n));
    }

    ThisType normalized() const { return ThisType(*this).normalize(); }

    virtual const MatrixRef<DType> T() const;
    virtual MatrixRef<DType> T();

    ThisType dot(const ThisType& rhs) const { return matmul(rhs); }

    template<typename RhsDType>
    auto matmul(const Matrix<RhsDType>& rhs) const -> Matrix<decltype(DType() * RhsDType())>
    {
        using ReturnType = Matrix<decltype(DType() * RhsDType())>;
        const ThisType& lhs(*this);
        if(lhs.shape(1) != rhs.shape(0))
        {
            // std::cout << "lhs.shape(1): " << lhs.shape(1) << ", rhs.shape(0): " << rhs.shape(0) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
        ReturnType ret = ReturnType::zeros({lhs.shape(0), rhs.shape(1)});
        ret.traverse([&](size_t i, size_t j)
            {
                for(int k = 0; k < lhs.shape(1); k++)
                    ret(i,j) += lhs(i, k) * rhs(k, j);
            });
        return ret;
    }

    std::string str() const { return to_string(*this); }

    ThisType& setBlock(size_t i0, size_t j0, const ThisType& mat)
    {
        mat.traverse([&](size_t i, size_t j) {(*this)(i + i0, j + j0) = mat(i, j);});
        return *this;
    }
    virtual const MatrixRef<DType> operator () (const Block& s) const;
    virtual MatrixRef<DType> operator () (const Block& s);
    ThisType pow(size_t n) const
    {
        ThisType ret = ThisType::ones(shape());
        for(size_t i = 0; i < n; i++) ret *= (*this);
        return ret;
    }


protected:
    Shape shape_;
    std::vector<DType> data_;
    bool major_;
};

using Mat = Matrix<FloatType>;

template<typename LType, typename DType> Matrix<DType> operator + (LType lhs, const Matrix<DType>& rhs) { return rhs + lhs;}
template<typename LType, typename DType> Matrix<DType> operator - (LType lhs, const Matrix<DType>& rhs) { return -rhs + lhs;}
#if 0
template<typename LType, typename DType> Matrix<DType> operator * (LType lhs, const Matrix<DType>& rhs) { return rhs * lhs;}
#else
template<typename LType, typename DType> typename std::enable_if_t< std::is_arithmetic<LType>::value , Matrix<DType>> operator * (LType lhs, const Matrix<DType>& rhs) { return rhs * lhs;}
template<typename LType, typename DType>
typename std::enable_if_t<
    std::is_same<
        typename mxm::Hypercomplex<typename LType::EntryType, LType::size()>, LType
    >::value, Matrix<DType>
>
operator * (LType lhs, const Matrix<DType>& rhs) { return rhs * lhs;}
#endif

inline size_t indexConvert2D(size_t i, size_t j, bool major, size_t shape_i, size_t shape_j)
{
    if(major == Mat::ROW) return (i * shape_j + j);
    return (j * shape_i + i);
}

inline void updateOffset(
    Shape& abs_offset, const Shape& inc_offset, bool same_major)
{
    if(same_major)
    {
        abs_offset[0] += inc_offset[0];
        abs_offset[1] += inc_offset[1];
    }else
    {
        abs_offset[0] += inc_offset[1];
        abs_offset[1] += inc_offset[0];
    }
}
#if 0
template <typename DType>
DType Matrix<DType>::norm(uint8_t p) const
{
    if(2 == p)
    {
        DType sum2(0);
        Matrix<DType> mat_2((*this)*(*this));
        traverse([&](size_t i, size_t j){sum2 += mat_2(i,j);});
        return sqrt(sum2);
    }
    if(1 == p)
    {
        DType max_abs_sum(0);
        for(int j = 0; j < shape(1); j++)
        {
            DType abs_sum(0);
            for(int i = 0; i < shape(0); i++) abs_sum += fabs((*this)(i,j));
            max_abs_sum = std::max(max_abs_sum, abs_sum);
        }
        return max_abs_sum;
    }

    throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    return 0.;
}
#endif

template<typename DType>
inline typename NormTraits<DType>::type norm(const Matrix<DType>& in)
{
    DType sum2(0);
    in.traverse([&](size_t i, size_t j){sum2 += in(i,j) * in(i,j);});
    return sqrt(sum2);
}

template<typename DType, unsigned int N>
DType norm(const Matrix<Hypercomplex<DType, N>>& in)
{
    DType sum2(0);
    in.traverse([&](size_t i, size_t j){ DType n = in(i,j).norm(); sum2 += n*n; });
    return sqrt(sum2);
}

template <typename DType>
typename NormTraits<DType>::type Matrix<DType>::norm(uint8_t p) const
{
    return mxm::norm(*this);
}

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, Matrix<DType>>::type
conj(const Matrix<DType>& in)
{
    return in;
}

template<typename DType, unsigned int N>
Matrix<Hypercomplex<DType, N>> conj(const Matrix<Hypercomplex<DType, N>>& in)
{
    Matrix<Hypercomplex<DType, N>> ret(in);
    ret.traverse([&](auto i, auto j){ret(i,j) = ret(i,j).conj();});
    return ret;
}

template<typename DType>
DType sum(const Matrix<DType>& src)
{
    DType ret(0);
    src.traverse([&](auto i, auto j){
        ret += src(i,j);
    });
    return ret;
}

} // namespace mxm

#ifdef MXM_HEADER_ONLY
// #include "linalg_mat_inl.h"
#endif

#endif // _LINALG_MAT_H
