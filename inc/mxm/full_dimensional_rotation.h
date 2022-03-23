#if !defined(_FULL_DIMENSIONAL_ROTATION_H_)
#define _FULL_DIMENSIONAL_ROTATION_H_

#include "linalg.h"
#include "rotation_base.h"
#include <iostream>


namespace mxm
{

template<typename DType, size_t DIM=3>
class Rotation
{
public:
    using ThisType = Rotation<DType, DIM>;

    Rotation(const Matrix<DType>& R): matrix_(R)
    {
        assert(matrix_.square());
        assert(matrix_.shape(0) == DIM);
    }

    template<typename RhsDType=DType>
    Rotation(const Rotation<RhsDType, DIM>& rhs)
    {
        matrix_ = rhs.asMatrix();
    }

    void operator = (const ThisType& rhs)
    {
        matrix_ = rhs.asMatrix();
    }

    template<typename RhsDType>
    void operator = (const Rotation<RhsDType, DIM>& rhs)
    {
        matrix_ = rhs.asMatrix();
    }

    // Rotation(size_t dim): matrix_(Matrix<DType>::identity(dim))
    // {
    // }

    Rotation(): matrix_(){}

    ThisType operator*(const ThisType& rhs) const { return ThisType(matrix_.matmul(rhs.matrix_)); }

    template<class DeriveType>
    Matrix<typename Traits<DeriveType>::EntryType>
    apply(const MatrixBase<DeriveType>& vector) const
    { return matrix_.matmul(vector); }
    // Matrix<DType> apply(const Matrix<DType>& vector) const { return matrix_.matmul(vector); }

    Matrix<DType> asMatrix() const { return matrix_; }

    ThisType inv() const { return ThisType(matrix_.T()); }

    constexpr size_t dim() const { return DIM; }

    template<class DeriveType> static
    std::enable_if_t<std::is_same<typename Traits<DeriveType>::EntryType, DType>::value, ThisType>
    fromMatrix(const MatrixBase<DeriveType>& R) { return ThisType(R); }
    static auto fromAngle(DType angle) { return Rotation<DType, 2>(mxm::rodrigues2D(angle)); }
    static auto fromAxisAngle(const Vector<DType>& axis, DType angle) { return Rotation<DType, 3>(mxm::rodrigues3D(axis, angle)); }
    static auto fromQuaternion(const Quaternion<DType>& q) { return Rotation<DType, 3>(mxm::toSO3(q)); }

    static ThisType fromPlaneAngle(const Vector<DType>& u, const Vector<DType>& v, DType angle)
    {
        auto bivec = planeAngleToBivector(u, v, angle);
        return ThisType(bivectorToRotationMatrix(bivec[0], bivec[1]));
    }

    static ThisType identity() { return ThisType(Matrix<DType>::identity(DIM)); }

    DType norm() const { return 0. ; }



private:
    Matrix<DType> matrix_;
};

} // namespace mxm


#endif // _FULL_DIMENSIONAL_ROTATION_H_
