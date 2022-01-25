#if !defined(_LINALG_MAT_REF_H_)
#define _LINALG_MAT_REF_H_

#include <vector>
#include <functional>
#include <numeric>
#include <array>
#include <string>
#include <assert.h>
#include <math.h>
#include <sstream>
#include <iomanip>

#include "linalg_matrix_base.h"


namespace mxm
{

// forward declaration
template<typename DType> class Matrix;
template<typename DType> class MatrixRef;
class Block;

template<typename DType>
struct Traits<MatrixRef<DType>>
{
    using EntryType = DType;
    using ArithType = typename Traits<DType>::ArithType; // norm type
    using DerefType = Matrix<DType>;
    static constexpr bool referable = true;
    // static std::string to_string(const MatrixRef<DType>& mat, size_t prec=6)
    // {
    //     return Traits<MatrixBase<MatrixRef<DType>>>::to_string(mat, prec);
    // }
};
void updateOffset(Shape& abs_offset, const Shape& inc_offset, bool same_major);
// end forward declaration

template<typename DType>
class MatrixRef: public MatrixBase<MatrixRef<DType>>
{
private:
    Shape shape_;
    Shape abs_offset_;
    Matrix<DType>* owner_;
    bool major_;

public:
    using ThisType = MatrixRef<DType>;
    using EntryType = typename Traits<ThisType>::EntryType;
    using ArithType = typename Traits<ThisType>::ArithType;

    // Constructors:
    // 01: Refer to standard Matrix type
    MatrixRef(Matrix<DType>& parent, const Shape& inc_offset, const Shape& shape, bool transpose);
    // 02: Refer to MatrixRef
    MatrixRef(MatrixRef<DType>& parent, const Shape& inc_offset, const Shape& shape, bool transpose);
    // 03: copy constructor
    MatrixRef(const ThisType& rhs);

    // assignment operators
    // 04 copy assignment
    template<typename DeriveType>  void operator = (const MatrixBase<DeriveType>& rhs);

    // 05 copy self type
    void operator = (const ThisType& rhs);

    //
    // static polymorphism
    const DType & operator () (size_t i, size_t j) const;
    DType & operator () (size_t i, size_t j);
    size_t shape(size_t ax) const {return shape_[ax];}
    const Shape& shape() const {return shape_;}

    MatrixRef<DType> operator () (const Block& b) { return MatrixBase<ThisType>::operator()(b); };
    const MatrixRef<DType> operator () (const Block& b) const { return MatrixBase<ThisType>::operator()(b); };

    // getters and setters
    const bool majorAxis() const {return major_;}
    const Matrix<DType>* owner() const { return owner_; }
};

//
// implementations
//

// 01:
template<typename DType>
MatrixRef<DType>::MatrixRef(Matrix<DType>& parent, const Shape& inc_offset, const Shape& shape, bool transpose)
:shape_(shape), abs_offset_{0,0}, owner_(&parent), major_(transpose ? !parent.majorAxis() : parent.majorAxis())
{
    updateOffset(abs_offset_, inc_offset, major_ == owner_->majorAxis());
}

// 02:
template<typename DType>
MatrixRef<DType>::MatrixRef(MatrixRef<DType>& parent, const Shape& inc_offset, const Shape& shape, bool transpose)
:shape_(shape), abs_offset_(parent.abs_offset_), owner_(parent.owner_), major_(transpose ? !parent.major_ : parent.major_)
{
    updateOffset(abs_offset_, inc_offset, major_ == owner_->majorAxis());
}

// 03: copy constructor
template<typename DType>
MatrixRef<DType>::MatrixRef(const ThisType& rhs)
: shape_(rhs.shape_), abs_offset_(rhs.abs_offset_), owner_(rhs.owner_), major_(rhs.major_)
{}

// 04: copy assignment
template<typename DType>
template<typename DeriveType>
void MatrixRef<DType>::operator = (const MatrixBase<DeriveType>& rhs_in)
{
    auto & rhs = reinterpret_cast<const DeriveType&>(rhs_in);
    assert(rhs.shape() == shape());

    this->traverse([&](auto i, auto j){ (*this)(i,j) = rhs(i,j); });

}

// 05: copy assignment
template<typename DType>
void MatrixRef<DType>::operator = (const MatrixRef<DType>& rhs)
{
    assert(rhs.shape() == shape());
    this->traverse([&](auto i, auto j){ (*this)(i,j) = rhs(i,j); });
}

template<typename DType>
const DType & MatrixRef<DType>::operator () (size_t i, size_t j) const
{
    Shape ref_ij(abs_offset_);
    updateOffset(ref_ij, {i,j}, owner_->majorAxis() == major_);
    return (*owner_)(ref_ij[0], ref_ij[1]);
}

template<typename DType>
DType & MatrixRef<DType>::operator () (size_t i, size_t j)
{
    return const_cast<DType &>(static_cast<const MatrixRef<DType>&>(*this)(i,j));
}

using MatRef = MatrixRef<float>;

} // namespace mxm



#endif // _LINALG_MAT_REF_H_
