#if !defined(_LINALG_MATRIX_BASE_H_)
#define _LINALG_MATRIX_BASE_H_

#include <vector>
#include <functional>
#include <numeric>
#include <array>
#include <string>
#include <assert.h>
#include <math.h>
#include <string>

namespace mxm
{
using Shape = std::array<size_t, 2>;
static const bool ROW = 0;
static const bool COL = 1;

template <typename DeriveType, typename=void>
struct Traits{};
template <typename DType> class Matrix;
template<typename DType> class MatrixRef;
class Block;
std::array<std::array<size_t, 2>, 2> deduct(const Block& b, const Shape& mat);

template<typename DType, typename=void>
typename Traits<DType>::ArithType norm(const DType& m);

template<typename DType>
struct Traits<DType, std::enable_if_t<std::is_arithmetic<DType>::value, void>>
{
    using ArithType = DType;
    static DType identity() { return DType(1); }
    static DType zero() { return DType(0); }
    static DType inv(DType val) { return identity() / val; }
    // static ArithType norm(DType val) { return abs(val); }
};


template <typename EntryType, typename ScalarType>
struct Operatable
{
    const static bool value = std::is_arithmetic<ScalarType>::value || std::is_same<ScalarType, EntryType>::value;
};

template<typename LhsType, typename RhsType, typename ReturnType_=decltype(LhsType()*RhsType())>
struct multiplies: public std::binary_function<LhsType,RhsType,ReturnType_>
{
    using ReturnType = ReturnType_;
    ReturnType operator()(const LhsType& lhs, const RhsType& rhs) { return lhs * rhs; }
};

template<typename LhsType, typename RhsType, typename ReturnType_=decltype(LhsType()*RhsType())>
struct plus: public std::binary_function<LhsType,RhsType,ReturnType_>
{
    using ReturnType = ReturnType_;
    ReturnType operator()(const LhsType& lhs, const RhsType& rhs) { return lhs + rhs; }
};

template<typename LhsType, typename RhsType, typename ReturnType_=decltype(LhsType()*RhsType())>
struct minus: public std::binary_function<LhsType,RhsType,ReturnType_>
{
    using ReturnType = ReturnType_;
    ReturnType operator()(const LhsType& lhs, const RhsType& rhs) { return lhs - rhs; }
};

//
// Curiously Recurring Template Pattern
//
template<typename DeriveType>
class MatrixBase
{
public:
    using ThisType = MatrixBase<DeriveType>;
    using EntryType = typename Traits<DeriveType>::EntryType;
    using ArithType = typename Traits<DeriveType>::ArithType;

    void traverse(std::function< void(size_t, size_t)> f) const
    {
        for(size_t i = 0; i < reinterpret_cast<const DeriveType*>(this)->shape(0); i++)
            for(size_t j = 0; j < reinterpret_cast<const DeriveType*>(this)->shape(1); j++)
                f(i,j);
    }

    EntryType& operator () (size_t i, size_t j) { return reinterpret_cast<DeriveType&>(*this)(i,j); }
    const EntryType& operator () (size_t i, size_t j) const { return reinterpret_cast<const DeriveType&>(*this)(i,j); }
    const Shape& shape() const { return reinterpret_cast<const DeriveType&>(*this).shape(); }
    size_t shape(size_t ax) const { return reinterpret_cast<const DeriveType&>(*this).shape(ax); }

    //
    // arithmetic operators
    //
#if 0
    template<typename RhsType, typename Op> std::enable_if_t<Operatable<EntryType, RhsType>::value , DeriveType>&
    inplaceOp(RhsType scalar, Op f)
    {
        auto& self = reinterpret_cast<DeriveType&>(*this);
        traverse([&](size_t i, size_t j){self(i,j) = f(self(i,j), scalar);});
        return self;
    }

    template<typename RhsType> DeriveType& operator *= (RhsType scalar) { return inplaceOp(scalar, multiplies<EntryType, RhsType>()); }
    template<typename RhsType> DeriveType& operator += (RhsType scalar) { return inplaceOp(scalar, plus<EntryType, RhsType>()); }
    template<typename RhsType> DeriveType& operator -= (RhsType scalar) { return inplaceOp(scalar, minus<EntryType, RhsType>()); }
#endif

    //  Compound Assignment Operator with Scalar
    template<typename RhsType> std::enable_if_t<Operatable<EntryType, RhsType>::value , DeriveType>&
    operator *= (RhsType scalar)  { auto& self = reinterpret_cast<DeriveType&>(*this); traverse([&](size_t i, size_t j){self(i,j) *= scalar;}); return self;}

    template<typename RhsType> std::enable_if_t<Operatable<EntryType, RhsType>::value , DeriveType>&
    operator += (RhsType scalar)  { auto& self = reinterpret_cast<DeriveType&>(*this); traverse([&](size_t i, size_t j){self(i,j) += scalar;}); return self;}

    template<typename RhsType> std::enable_if_t<Operatable<EntryType, RhsType>::value , DeriveType>&
    operator -= (RhsType scalar)  { auto& self = reinterpret_cast<DeriveType&>(*this); traverse([&](size_t i, size_t j){self(i,j) -= scalar;}); return self;}


    // Binary Operator
    template<typename RhsType> std::enable_if_t<Operatable<EntryType, RhsType>::value , Matrix<EntryType>>
    operator * (RhsType scalar) const { return Matrix<EntryType>(reinterpret_cast<const DeriveType&>(*this)) *= scalar; }

    template<typename RhsType> std::enable_if_t<Operatable<EntryType, RhsType>::value , Matrix<EntryType>>
    operator + (RhsType scalar) const { return Matrix<EntryType>(reinterpret_cast<const DeriveType&>(*this)) += scalar; }

    template<typename RhsType> std::enable_if_t<Operatable<EntryType, RhsType>::value , Matrix<EntryType>>
    operator - (RhsType scalar) const { return Matrix<EntryType>(reinterpret_cast<const DeriveType&>(*this)) -= scalar; }


    template<typename T> DeriveType & operator *= (const MatrixBase<T>& rhs)  { auto& self = reinterpret_cast<DeriveType&>(*this); traverse([&](size_t i, size_t j){self(i,j) *= reinterpret_cast<const T &>(rhs)(i,j);}); return self;}
    template<typename T> DeriveType & operator += (const MatrixBase<T>& rhs)  { auto& self = reinterpret_cast<DeriveType&>(*this); traverse([&](size_t i, size_t j){self(i,j) += reinterpret_cast<const T &>(rhs)(i,j);}); return self;}
    template<typename T> DeriveType & operator -= (const MatrixBase<T>& rhs)  { auto& self = reinterpret_cast<DeriveType&>(*this); traverse([&](size_t i, size_t j){self(i,j) -= reinterpret_cast<const T &>(rhs)(i,j);}); return self;}

    template<typename T, typename OpType>
    auto binaryOp(const MatrixBase<T>& rhs, OpType f) const
    {
        auto& self = reinterpret_cast<const DeriveType&>(*this);
        auto& rhs_d = reinterpret_cast<const T&>(rhs);
        Matrix<typename OpType::ReturnType> ret(self.shape());
        ret.traverse([&](size_t i, size_t j){ret(i,j) = f(self(i,j), rhs_d(i,j));});
        return ret;
    }
#if 1
    template<typename T> auto operator * (const MatrixBase<T>& rhs) const { return binaryOp(rhs, multiplies<EntryType, typename Traits<T>::EntryType>()); }
    template<typename T> auto operator + (const MatrixBase<T>& rhs) const { return binaryOp(rhs, plus<EntryType, typename Traits<T>::EntryType>()); }
    template<typename T> auto operator - (const MatrixBase<T>& rhs) const { return binaryOp(rhs, minus<EntryType, typename Traits<T>::EntryType>()); }
#else
    template<typename T> auto operator * (const MatrixBase<T>& rhs)
    {
        auto& self = reinterpret_cast<DeriveType&>(*this);
        Matrix<decltype(typename Traits<DeriveType>::EntryType() * typename Traits<T>::EntryType())> ret(self.shape());
        ret.traverse([&](size_t i, size_t j){ret(i,j) = self(i,j) * reinterpret_cast<const T&>(rhs)(i,j);});
        return ret;
    }
    template<typename T> auto operator + (const MatrixBase<T>& rhs)
    {
        auto& self = reinterpret_cast<DeriveType&>(*this);
        Matrix<decltype(typename Traits<DeriveType>::EntryType() + typename Traits<T>::EntryType())> ret(self.shape());
        ret.traverse([&](size_t i, size_t j){ret(i,j) = self(i,j) + reinterpret_cast<const T&>(rhs)(i,j);});
        return ret;
    }
    template<typename T> auto operator - (const MatrixBase<T>& rhs)
    {
        auto& self = reinterpret_cast<DeriveType&>(*this);
        Matrix<decltype(typename Traits<DeriveType>::EntryType() - typename Traits<T>::EntryType())> ret(self.shape());
        ret.traverse([&](size_t i, size_t j){ret(i,j) = self(i,j) - reinterpret_cast<const T&>(rhs)(i,j);});
        return ret;
    }
#endif

    DeriveType operator -() const        { return reinterpret_cast<DeriveType&>(*this) *= -1;}

    bool operator == (const DeriveType& rhs)
    {
        auto& self = reinterpret_cast<DeriveType&>(*this);
        if(rhs.shape() != self.shape()) return false;

        for(size_t i = 0; i < self.shape(0); i++)
            for(size_t j = 0; j < self.shape(1); j++)
                if(self(i,j) != rhs(i,j)) return false;
        return true;
    }
    bool operator != (const DeriveType& rhs) { return ! (*this == rhs); }

    // end of operator overloading
    //

    template<typename T>
    auto matmul(const MatrixBase<T>& rhs) const
    {
        auto& self = reinterpret_cast<const DeriveType&>(*this);
        const T& d_rhs = reinterpret_cast<const T&>(rhs);


        assert(self.shape(1) == d_rhs.shape(0) && "matrix multiplication shape mismatch");

        Matrix<decltype(typename Traits<DeriveType>::EntryType() * typename Traits<T>::EntryType())> ret({self.shape(0), d_rhs.shape(1) });
        ret.traverse([&](size_t i, size_t j)
            {
                ret(i,j) = 0;
                for(int k = 0; k < self.shape(1); k++)
                    ret(i,j) += self(i, k) * d_rhs(k, j);
            });
        return ret;
    }

    template<typename T>
    DeriveType& setBlock(size_t i0, size_t j0, const MatrixBase<T>& mat)
    {
        mat.traverse([&](size_t i, size_t j) {(*this)(i + i0, j + j0) = reinterpret_cast<const T&>(mat)(i, j);});
        return reinterpret_cast<DeriveType&>(*this);
    }

    //Sub Matrix by Copy
    template<bool Copy=false> std::enable_if_t<Copy, Matrix<EntryType>> T() const;

    // Sub Matrix by Ref
    template<bool Copy=false> std::enable_if_t<!Copy && Traits<DeriveType>::referable, const MatrixRef<EntryType>> T() const;
    template<bool Copy=false> std::enable_if_t<!Copy && Traits<DeriveType>::referable, MatrixRef<EntryType>> T();
    std::enable_if_t<Traits<DeriveType>::referable , MatrixRef<EntryType>> operator () (const Block& b);
    std::enable_if_t<Traits<DeriveType>::referable , const MatrixRef<EntryType>> operator () (const Block& b) const;

#if 1
    // static methods
    static DeriveType zeros(const Shape& shape) { return DeriveType(shape) *= Traits<EntryType>::zero() ; }
    static DeriveType ones(const Shape& shape) { return zeros(shape) += Traits<EntryType>::identity(); };
    static DeriveType identity(size_t n)
    {
        DeriveType ret = zeros({n,n});
        for(size_t i = 0; i < n; i++)
            ret(i,i) = Traits<EntryType>::identity();
        return ret;
    }
#endif

    bool square() const {
        auto& self = reinterpret_cast<const DeriveType&>(*this);
        return self.shape(0) == self.shape(1);
    }
    EntryType trace() const
    {
        auto& self = reinterpret_cast<const DeriveType&>(*this);
        assert(square() && "Only square matrix has trace");
        EntryType ret = Traits<EntryType>::zero();
        for(size_t i = 0;i < self.shape(0); i++) ret += self(i,i);
        return ret;
    }

    ArithType norm() const;
    DeriveType& normalize();
    Matrix<EntryType> normalized() const;
    Matrix<EntryType> inv() const;

private:
};

template<typename ScalarType, typename DeriveType>
std::enable_if_t<Operatable<typename Traits<DeriveType>::EntryType, ScalarType>::value , typename Traits<DeriveType>::DerefType>
operator * (ScalarType scalar, const MatrixBase<DeriveType>& rhs) { return rhs * scalar; }

template<typename ScalarType, typename DeriveType>
std::enable_if_t<Operatable<typename Traits<DeriveType>::EntryType, ScalarType>::value , typename Traits<DeriveType>::DerefType>
operator + (ScalarType scalar, const MatrixBase<DeriveType>& rhs) { return rhs + scalar; }

template<typename ScalarType, typename DeriveType>
std::enable_if_t<Operatable<typename Traits<DeriveType>::EntryType, ScalarType>::value , typename Traits<DeriveType>::DerefType>
operator - (ScalarType scalar, const MatrixBase<DeriveType>& rhs) { return rhs - scalar; }

#if 0
template <typename T, typename=void>
struct ArithTraits{};

template <typename T>
struct ArithTraits<T, std::enable_if_t<std::is_arithmetic<T>::value, void>>{using type = T;};
#endif

//
// implementations
//

template<typename DeriveType>
template<bool Copy>
std::enable_if_t<Copy, Matrix<typename Traits<DeriveType>::EntryType>> MatrixBase<DeriveType>::T() const
{
    auto & self = reinterpret_cast<const DeriveType&>(*this);

    Matrix<EntryType> ret({self.shape(1), self.shape(0)});
    ret.traverse([&](auto i, auto j){ ret(i,j) = self(i,j); });
    return ret;
}

template<typename DeriveType>
template<bool Copy>
std::enable_if_t<!Copy && Traits<DeriveType>::referable, const MatrixRef<typename Traits<DeriveType>::EntryType>>
MatrixBase<DeriveType>::T() const
{
    auto & self = reinterpret_cast<const DeriveType&>(*this);
    return MatrixRef<EntryType>(const_cast<DeriveType&>(self), {0,0}, {self.shape(1), self.shape(0)}, true);
}

template<typename DeriveType>
template<bool Copy>
std::enable_if_t<!Copy && Traits<DeriveType>::referable, MatrixRef<typename Traits<DeriveType>::EntryType>>
MatrixBase<DeriveType>::T()
{
    auto & self = reinterpret_cast<DeriveType&>(*this);
    return MatrixRef<EntryType>(self, {0,0}, {self.shape(1), self.shape(0)}, true);
}

template<typename DeriveType>
std::enable_if_t<Traits<DeriveType>::referable , const MatrixRef<typename Traits<DeriveType>::EntryType>>
MatrixBase<DeriveType>::operator () (const Block& b) const
{
    auto & self = reinterpret_cast<const DeriveType&>(*this);
    auto offset_shape = deduct(b, self.shape());
    return MatrixRef<EntryType>(const_cast<DeriveType&>(self), offset_shape[0], offset_shape[1], false);
}

template<typename DeriveType>
std::enable_if_t<Traits<DeriveType>::referable , MatrixRef<typename Traits<DeriveType>::EntryType>>
MatrixBase<DeriveType>::operator () (const Block& b)
{
    auto & self = reinterpret_cast<DeriveType&>(*this);
    auto offset_shape = deduct(b, self.shape());
    return MatrixRef<EntryType>(self, offset_shape[0], offset_shape[1], false);
}

#if 0

template<typename DType, typename=void>
std::string to_string(const DType& m, size_t prec=6);

template<typename DeriveType, typename=void>
std::string to_string(const MatrixBase<DeriveType>& mat_in, size_t prec=6)
{
    auto & mat = reinterpret_cast<const DeriveType&>(mat_in);
    std::string ret;
    mat.traverse([&](size_t i, size_t j){
        ret += (mxm::to_string(mat(i,j), prec) + (j == mat.shape(1) - 1 ? "\n" : " ")); });
    return ret;
}

// template<typename DType, typename=std::enable_if_t<std::is_arithmetic<DType>::value>>
template<typename DType, typename=std::enable_if_t<std::is_same<float, DType>::value>>
std::string
to_string(const DType& v, size_t prec)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(prec) << v;
    return stream.str();
}

#endif

template<typename DeriveType>
struct Traits<MatrixBase<DeriveType>>
{
    using EntryType = typename Traits<DeriveType>::EntryType ;
    using ArithType = typename Traits<EntryType>::ArithType; // norm type
    static constexpr bool referable = true;
};
// Inversion

template<typename DeriveType, typename=void>
Matrix<typename Traits<DeriveType>::EntryType>
inv(const MatrixBase<DeriveType>& mat);

template<typename DeriveType>
Matrix< typename Traits<DeriveType>::EntryType> MatrixBase<DeriveType>::inv() const
{
    return mxm::inv(*this);
}

// Norm and normalization interfaces
template<typename DeriveType, typename=void>
typename Traits<DeriveType>::ArithType
norm(const MatrixBase<DeriveType>& mat)
{
    using EntryType = typename Traits<DeriveType>::EntryType;
    using ArithType = typename Traits<DeriveType>::ArithType;
    auto& self = reinterpret_cast<const DeriveType&>(mat);
    typename Traits<DeriveType>::ArithType sum2(0);
    self.traverse([&](size_t i, size_t j){
        auto n = norm(self(i,j)); sum2 += n*n;
    });
    return sqrt(sum2);
}

template<typename DeriveType>
typename Traits<DeriveType>::ArithType MatrixBase<DeriveType>::norm() const
{
    return mxm::norm(*this);
}

template<typename DeriveType>
DeriveType& MatrixBase<DeriveType>::normalize()
{
    auto & self = reinterpret_cast<DeriveType&>(*this);
    self *= Traits<ArithType>::inv(mxm::norm(*this));
    return self;
}

template<typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType> MatrixBase<DeriveType>::normalized() const
{
    using EntryType = typename Traits<DeriveType>::EntryType;
    Matrix<EntryType> ret(*this);
    return ret.normalize();
}

// for float point type

template<typename DType,
    typename=std::enable_if_t<std::is_floating_point<DType>::value, void>>
typename Traits<DType>::ArithType norm(const DType& val){ return std::abs(val); }
template<typename DType,
    typename=std::enable_if_t<std::is_floating_point<DType>::value, void>>
typename Traits<DType>::ArithType inv(const DType& val){ return DType(1)/val; }

template<typename DeriveType, typename=std::enable_if_t<std::is_arithmetic<typename Traits<DeriveType>::EntryType>::value, void>>
Matrix<typename Traits<DeriveType>::EntryType>
conj(const MatrixBase<DeriveType>& in)
{
    return in;
}
#if 0
template<typename ArithType, unsigned int N>
class Hypercomplex;

template<typename DeriveType,
    typename=std::enable_if_t<
        std::is_same<
            Hypercomplex<
                typename Traits<DeriveType>::ArithType,
                Traits<DeriveType>::EntryType::size()
            >, typename Traits<DeriveType>::EntryType
        >::value, void
    >
>
Matrix<Hypercomplex<typename Traits<DeriveType>::ArithType, Traits<DeriveType>::EntryType::size()>>
conj(const MatrixBase<DeriveType>& in)
{
    using ArithType = typename Traits<DeriveType>::ArithType;
    const size_t N = Traits<DeriveType>::EntryType::size();
    Matrix<Hypercomplex<ArithType, N>> ret(in);
    ret.traverse([&](auto i, auto j){ret(i,j) = ret(i,j).conj();});
    return ret;
}
#endif
} // namespace mxm



#endif // _LINALG_MATRIX_BASE_H_
