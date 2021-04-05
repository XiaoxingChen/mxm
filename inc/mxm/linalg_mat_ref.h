#if !defined(_LINALG_MAT_REF_H)
#define _LINALG_MAT_REF_H

#include "linalg_mat.h"


namespace mxm
{

template<typename DType>
class VectorRef;

template<typename DType>
class MatrixRef: public Matrix<DType>
{
public:
    using ThisType = MatrixRef;
    using BaseType = Matrix<DType>;

    MatrixRef(BaseType& parent, const Shape& inc_offset, const Shape& shape, bool transpose)
        :BaseType({0,0}, {}, transpose ? !parent.majorAxis() : parent.majorAxis()),
        owner_(parent.owner()), abs_offset_(parent.absOffset())
    {
        updateOffset(abs_offset_, inc_offset, BaseType::majorAxis() == owner().majorAxis());
        BaseType::shape_ = shape;
    }

    MatrixRef(const ThisType& rhs)
        :BaseType({0,0}, {}, rhs.majorAxis()),
        abs_offset_(rhs.absOffset()), owner_(const_cast<BaseType&>(rhs.owner()))
    {
        BaseType::shape_ = rhs.shape();
    }


    virtual BaseType& owner() override { return owner_; }
    virtual const BaseType& owner() const override { return owner_; }
    virtual const Shape& absOffset() const override { return abs_offset_; }
    virtual const DType& operator () (size_t i, size_t j) const override
    {
        Shape ref_ij = absOffset();
        updateOffset(ref_ij, {i,j}, owner().majorAxis() == BaseType::majorAxis());
        return owner_(ref_ij[0], ref_ij[1]);
    }

    virtual DType& operator () (size_t i, size_t j) override { return const_cast<DType&>(static_cast<const ThisType&>(*this)(i,j)); }
    virtual const MatrixRef<DType> operator () (const Block& s) const override
    {
        return Matrix<DType>::operator () (s);
    }

    virtual MatrixRef<DType> operator () (const Block& s) override
    {
        return static_cast<const ThisType&>(*this)(s);
    }

    //
    // All overloaded operators except assignment (operator=)
    // are inherited by derived classes.
    virtual void operator = (const BaseType& rhs) override
    {
        if(BaseType::shape() != rhs.shape())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        Matrix<DType>::traverse([&](size_t i, size_t j) { (*this)(i,j) = rhs(i,j); });
    }

    void operator = (const ThisType& rhs) { (*this) = static_cast<const BaseType&>(rhs); }

    const VectorRef<DType> asVector() const;
    VectorRef<DType> asVector();

protected:
    BaseType& owner_;
    Shape abs_offset_;

};

template <typename DType>
const MatrixRef<DType> Matrix<DType>::T() const
{
    return MatrixRef<DType>(
        const_cast<ThisType&>(*this), {0,0},
        {shape(1), shape(0)}, true);
}

template <typename DType>
MatrixRef<DType> Matrix<DType>::T()
{
    return static_cast<const Matrix<DType>&>(*this).T();
}

template <typename DType>
const MatrixRef<DType> Matrix<DType>::operator () (const Block& s) const
{
    auto i01j01 = s.getBlock(*this);
    return MatrixRef<DType>(
        const_cast<ThisType&>(*this),
        {i01j01[0],i01j01[2]},
        {i01j01[1] - i01j01[0], i01j01[3] - i01j01[2]}, false);
}

template <typename DType>
MatrixRef<DType> Matrix<DType>::operator () (const Block& s)
{
    return static_cast<const Matrix<DType>&>(*this)(s);
}

using MatRef = MatrixRef<FloatType>;


} // namespace mxm

#endif // _LINALG_MAT_REF_H
