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

#include "linalg_matrix_base.h"

namespace mxm
{

// forward declaration
size_t index1D(size_t i, size_t j, bool major, const Shape& shape);
template<typename DType> class Matrix;
template<typename DType> class MatrixRef;
class Block;
std::array<std::array<size_t, 2>, 2> deduct(const Block& b, const Shape& mat);

class AutoShape;
AutoShape fixRow(size_t n);
AutoShape fixCol(size_t n);
AutoShape square();

template<typename DType>
struct Traits<Matrix<DType>>
{
    using EntryType = DType;
    using ArithType = typename Traits<DType>::ArithType; // norm type
    using DerefType = Matrix<DType>;
    static constexpr bool referable = true;
    // static std::string to_string(const Matrix<DType>& mat, size_t prec=6)
    // {
    //     return Traits<MatrixBase<Matrix<DType>>>::to_string(mat, prec);
    // }
};
// end forward declaration

template<typename DType>
class Matrix: public MatrixBase<Matrix<DType>>
{
private:
    Shape shape_;
    std::vector<DType> data_;
    bool major_;

public:
    using ThisType = Matrix<DType>;
    using EntryType = typename Traits<ThisType>::EntryType;
    using ArithType = typename Traits<ThisType>::ArithType;

    ~Matrix(){}

    // Constructors:
    // Constructor 00: default
    Matrix();
    // Constructor 01: move initialize
    Matrix(const AutoShape& _shape, std::vector<DType>&& data, bool major=ROW);
    // Constructor 02:
    Matrix(const AutoShape& _shape, const std::vector<DType>& data={}, bool major=ROW);
    // Constructor 03: standard copy constructor
    Matrix(const ThisType& rhs);
    // Constructor 04: standard move constructor
    Matrix(ThisType&& rhs);
    // Constructor 05: type cast constructor
    template<typename DeriveType> Matrix(const MatrixBase<DeriveType>& rhs);
    // Constructor 06: construct from MatrixRef
    Matrix(const MatrixRef<DType>& rhs);

    // assignment operators
    // copy assignment
    void operator = (const ThisType& rhs);
    // move assignment
    void operator = (ThisType&& rhs);


    //
    // static polymorphism
    const DType & operator () (size_t i, size_t j) const { return data_.at(index1D(i,j,major_, shape_)); }
    DType & operator () (size_t i, size_t j) { return data_.at(index1D(i,j,major_, shape_)); }
    size_t shape(size_t ax) const {return shape_[ax];}
    const Shape& shape() const {return shape_;}

    //
    const bool majorAxis() const {return major_;}

    MatrixRef<DType> operator () (const Block& b) { return MatrixBase<ThisType>::operator()(b); };
    const MatrixRef<DType> operator () (const Block& b) const { return MatrixBase<ThisType>::operator()(b); };

    const DType * data() const { return data_.data(); }
    void reshape(const AutoShape& shape, bool major);

};

class AutoShape
{
public:
    enum EnumDynamic{
        eRowDefined = 0,
        eColDefined = 1,
        eFullyDefined,
        eSquare,
        eUndefined};
#if 0
    AutoShape(const Shape& shape): shape_(shape), defined_{true, true} {}
    AutoShape(std::initializer_list<size_t> shape): shape_{*shape.begin(), *(shape.begin() + 1)}, defined_{true, true} {}
    AutoShape(size_t row, size_t col, bool row_def, bool col_def): shape_{row, col}, defined_{row_def, col_def} {}

#else
    AutoShape(const Shape& shape): shape_(shape), state_(eFullyDefined) {}
    AutoShape(std::initializer_list<size_t> shape): shape_{*shape.begin(), *(shape.begin() + 1)}, state_(eFullyDefined) { assert(2 == shape.size()); }
    AutoShape(const EnumDynamic& state, size_t row, size_t col): shape_{row, col}, state_(state) {}
#endif
    std::array<size_t, 2> deduct(size_t total_num) const;

private:
    std::array<size_t, 2> shape_;
    // std::array<bool, 2> defined_;
    EnumDynamic state_;
};

//
// implementations
//

// 00
template<typename DType>
Matrix<DType>::Matrix(): shape_({0,0}), data_(0), major_(ROW)
{

}

// 01
template<typename DType>
Matrix<DType>::Matrix(const AutoShape& _shape, std::vector<DType>&& data, bool major)
    :shape_(_shape.deduct(data.size())), data_(std::move(data)), major_(major) {}

// 02
template<typename DType>
Matrix<DType>::Matrix(const AutoShape& _shape, const std::vector<DType>& data, bool major)
    :shape_(_shape.deduct(data.size())), data_(shape_[0] * shape_[1]), major_(major)
{
    if(data.size() > 0)
        data_ = data;
}

// 03
template<typename DType>
Matrix<DType>::Matrix(const ThisType& rhs)
:shape_(rhs.shape_), data_(rhs.data_), major_(rhs.major_)
{}

// 04
template<typename DType>
Matrix<DType>::Matrix(ThisType&& rhs)
:shape_(rhs.shape_), data_(std::move(rhs.data_)), major_(rhs.major_)
{
    rhs.shape_ = {0,0};
}

// 05
template<typename DType>
template<typename DeriveType>
Matrix<DType>::Matrix(const MatrixBase<DeriveType>& rhs)
    :shape_(reinterpret_cast<const DeriveType&>(rhs).shape()), data_(shape_[0] * shape_[1]), major_(ROW)
{
    this->traverse([&](auto i, auto j){
        (*this)(i,j) = DType(reinterpret_cast<const DeriveType&>(rhs)(i,j));
    });
}

// 06
template<typename DType>
Matrix<DType>::Matrix(const MatrixRef<DType>& rhs)
{
    if(rhs.owner() == this) // self assignment
    {
        if(rhs.majorAxis() != major_ && rhs.shape() == Shape{shape_[1], shape_[0]})
        {//deal with `mat = mat.T();`
            major_ = ! major_;
        }else if(rhs.majorAxis() == major_ && rhs.shape() == shape_)
        {// deal with `mat = mat(Block({0,end()},{0,end()}));
            return ;
        }else // deal with general cases like `mat = mat(Block({1,end()}, {1,end()})).T();`
        {     // create temporary object to avoid memory overlap

            ThisType tmp(rhs);
            (*this) = std::move(tmp);
        }
    }else // rhs is unrelated with this
    {
        shape_ = rhs.shape();
        data_.resize(shape_[0] * shape_[1]);
        major_ = rhs.majorAxis();
        this->traverse([&](size_t i, size_t j){ (*this)(i,j) = rhs(i,j); });
    }
}

template<typename DType>
void Matrix<DType>::operator = (const Matrix<DType>& rhs)
{
    if(&rhs == this) return; // deal with `mat = mat;`
    shape_ = rhs.shape();
    data_ = rhs.data_;
    major_ = rhs.major_;
}

template<typename DType>
void Matrix<DType>::operator = (Matrix<DType>&& rhs)
{
    if(&rhs == this) return; // deal with `mat = mat;`
    shape_ = rhs.shape();
    data_.swap(rhs.data_);
    major_ = rhs.major_;
    rhs.shape_ = {0,0};
}
//
// end of ref operation

template<typename DType>
void Matrix<DType>::reshape(const AutoShape& shape, bool major)
{
    auto new_shape = shape.deduct(data_.size());
    assert(new_shape[0] * new_shape[1] == data_.size());
    shape_ = new_shape;
    major_ = major;
}

inline size_t index1D(size_t i, size_t j, bool major, const Shape& shape)
{
    if(major == ROW) return (i * shape[1] + j);
    return (j * shape[0] + i);
}

template<typename EntryType, typename=void>
typename Matrix<EntryType>::ArithType
norm(const Matrix<EntryType>& mat)
{
    using BaseType = MatrixBase<Matrix<EntryType>>;
    return norm(reinterpret_cast<const BaseType&>(mat));
}
#if 0
template<typename EntryType, typename=void>
std::string
to_string(const Matrix<EntryType>& mat, size_t prec=6)
{
    using BaseType = MatrixBase<Matrix<EntryType>>;
    return to_string(reinterpret_cast<const BaseType&>(mat), prec);
}
#endif

#if 0
template<typename EntryType, typename=void>
Matrix<EntryType>&
normalize(Matrix<EntryType>& mat)
{
    using ArithType = typename Matrix<EntryType>::ArithType;
    mat *= Traits<ArithType>::inv(mxm::norm(mat));
    return mat;
}

template<typename EntryType, typename=void>
Matrix<EntryType>
normalized(const Matrix<EntryType>& mat)
{
    Matrix<EntryType> ret(mat);
    return normalize(ret);
}
#endif


using Mat = Matrix<float>;
} // namespace mxm

#ifdef MXM_HEADER_ONLY
// #include "linalg_mat_inl.h"
#endif

#endif // _LINALG_MAT_H