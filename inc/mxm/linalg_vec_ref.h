#if !defined(_LINALG_VEC_REF_H)
#define _LINALG_VEC_REF_H

#include "mxm/linalg_vec.h"
#include "mxm/linalg_mat_block.h"


namespace mxm
{

template<typename DType>
class VectorRef: public Vector<DType>
{
public:
    using ThisType = VectorRef;
    using BaseType = Vector<DType>;
    using Base2Type = Matrix<DType>;

    // parent: can be either Matrix or Vector
    VectorRef(Base2Type& parent, const Shape& inc_offset, const Shape& shape, bool transpose)
        :BaseType(), owner_(parent.owner()), abs_offset_(parent.absOffset())
    {
        // if(BaseType::shape(0) != 1 && BaseType::shape(1) != 1)
        //     throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        updateOffset(abs_offset_, inc_offset, BaseType::majorAxis() == BaseType::owner().majorAxis());
        Base2Type::shape_ = shape;
        Base2Type::major_ = transpose ? !parent.owner().majorAxis() : parent.owner().majorAxis();
    }

    VectorRef(const ThisType& rhs)
        :BaseType(), owner_(const_cast<Base2Type&>(rhs.owner())), abs_offset_(rhs.absOffset())
    {
        Base2Type::shape_ = rhs.shape();
        Base2Type::major_ = rhs.majorAxis();
    }

    //
    // All overloaded operators except assignment (operator=)
    // are inherited by derived classes.
    void operator = (const BaseType& rhs)
    {
        Base2Type::operator = (rhs);
    }

    void operator = (const ThisType& rhs) { (*this) = static_cast<const BaseType&>(rhs); }

    virtual const DType& operator () (size_t i) const override
    {
        if(1 == BaseType::shape(0))
            return (*this)(0, i);
        return (*this)(i, 0);
    }

    virtual DType& operator () (size_t i) override
    {
        return const_cast<DType&>(static_cast<const ThisType&>(*this)(i));
    }

    virtual const DType& operator () (size_t i, size_t j) const override
    {
        Shape ref_ij = abs_offset_;
        updateOffset(ref_ij, {i,j}, BaseType::owner().majorAxis() == BaseType::majorAxis());
        return owner_(ref_ij[0], ref_ij[1]);
    }

    virtual DType& operator () (size_t i, size_t j) override { return const_cast<DType&>(static_cast<const ThisType&>(*this)(i,j)); }

    virtual const Shape& absOffset() const override { return abs_offset_; }
    virtual Base2Type& owner() override { return owner_; }
    virtual const Base2Type& owner() const override { return owner_; }

private:
    Base2Type& owner_;
    Shape abs_offset_;
};

template <typename DType>
const VectorRef<DType> MatrixRef<DType>::asVector() const
{
    if(BaseType::shape(0) != 1 && BaseType::shape(1) != 1)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    using ReturnType = VectorRef<DType>;
    return static_cast<const ReturnType>(ReturnType(
            const_cast<ThisType&>(*this),
            /*offset*/ {0, 0},
            /*shape*/ BaseType::shape(),
            /*transpose*/ false));
}

template <typename DType>
VectorRef<DType> MatrixRef<DType>::asVector()
{
    using ReturnType = VectorRef<DType>;
    return static_cast<const ThisType&>(*this).asVector();
}

using VecRef = VectorRef<FloatType>;

} // namespace mxm



#endif // _LINALG_VEC_REF_H
