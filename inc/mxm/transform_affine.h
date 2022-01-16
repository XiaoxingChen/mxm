#if !defined(__TRANSFORM_TRS_H__)
#define __TRANSFORM_TRS_H__

#include "linalg.h"
#include "full_dimensional_rotation.h"

namespace mxm
{

// U = K S
// U is the upper triangle matrix
// K is the diagonal matrix that indicates the scale
// S is the shear matrix
template<typename DType>
Vector<DType> scaleFromUpperTriangle(const Matrix<DType>& upper_tri)
{

    assert(upper_tri.square());
    size_t dim = upper_tri.shape(0);
    Vector<DType> scale(dim);
    auto norm2 = mxm::sum(upper_tri * upper_tri, 1);
    for(size_t i = 0; i < dim; i++)
    {
        scale(i) = sqrt(norm2(i, 0));
    }
    return scale;
}

template<typename DType>
Matrix<DType> simpleReflectMatrix_(size_t dim)
{
    auto mat = Matrix<DType>::identity(dim);
    Matrix<DType> tmp = mat(Col(0));
    mat(Col(0)) = mat(Col(1));
    mat(Col(1)) = tmp;
    return mat;
}

// References: https://math.stackexchange.com/questions/1120209/decomposition-of-4x4-or-larger-affine-transformation-matrix-to-individual-variab
//  T * R * S * Sh
template<typename DType, size_t DIM=3>
class AffineTransform
{
public:
    using ThisType = AffineTransform<DType, DIM>;

    AffineTransform(
        const Vector<DType>& trans,
        const Rotation<DType, DIM>& r,
        const Vector<DType>& scale,
        const Matrix<DType>& shear=Matrix<DType>::identity(DIM))
        :translation_(trans), rotation_(r), scale_(scale), shear_(shear)
    {
        assert(trans.size() == dim());
        assert(scale.size() == dim());
        assert(shear.square());
        assert(shear.shape(0) == dim());
    }

    AffineTransform():
        translation_(Vector<DType>::zeros(DIM)),
        rotation_(Rotation<DType, DIM>::identity()),
        scale_(Vector<DType>::ones(DIM)),
        shear_(Matrix<DType>::identity(DIM))
    {}

    ~AffineTransform(){};

    Matrix<DType> apply(const Matrix<DType>& vectors) const
    {
        // T * R * S * Sh
        auto ret = rotation_.apply( mxm::diagonalMatrix(scale_).matmul(shear_).matmul(vectors));
        ret += translation_;
        return ret;
    }

    constexpr size_t dim() const { return DIM; }

    ThisType inv() const {
        return ThisType::fromMatrix(mxm::inv(asMatrix()));
    }

    Matrix<DType> asMatrix() const
    {
        Matrix<DType> ret = Matrix<DType>::identity(dim() + 1);
        ret.setBlock(0,0, rotation_.asMatrix().matmul (mxm::diagonalMatrix(scale_)).matmul(shear_) );
        ret.setBlock(0, DIM, translation_);
        return ret;
    }

    // ThisType inv() const { return ThisType(rotation_.inv().apply(-translation_), rotation_.inv()); }

    ThisType operator*(const ThisType& rhs) const
    {
        return fromMatrix(asMatrix().matmul(rhs.asMatrix()));
    }

    static ThisType identity()
    {
        return ThisType();
    }

    static ThisType fromMatrix(const Matrix<DType>& mat)
    {
        //affine decomposition
        auto q_r = mxm::qr::decomposeByRotation(mat(Block({0,DIM}, {0, DIM})));
        Matrix<DType>& rot(q_r[0]);
        Matrix<DType>& upper_tri(q_r[1]);
        {
            // check reflection
            auto expect_eye = rot.matmul(rot.T());
            DType det(1);
            for(size_t i = 0; i < DIM; i++) det *= expect_eye(i,i);
            if(det < 0)
            {
                auto reflect = simpleReflectMatrix_<DType>(DIM);
                rot = rot.matmul(reflect);
                upper_tri = reflect.matmul(upper_tri);
            }
        }
        auto scale_vec = scaleFromUpperTriangle(upper_tri);
        Vector<DType> inv_scale_vec(scale_vec);
        for(size_t i = 0; i < DIM; i++) inv_scale_vec(i) = DType(1) / scale_vec(i);
        Matrix<DType> shear = mxm::diagonalMatrix(inv_scale_vec).matmul(upper_tri);
        auto tra = Vector<DType>(mat(Block({0, DIM}, {DIM, DIM + 1})));

        return ThisType(tra, Rotation<DType, DIM>::fromMatrix(rot), scale_vec, shear);
    }

    const Rotation<DType, DIM>& rotation() const {return rotation_;}
    Rotation<DType, DIM>& rotation() {return rotation_;}
    const Vector<DType>& translation() const {return translation_;}
    Vector<DType>& translation() {return translation_;}
    const Vector<DType>& scale() const {return scale_;}
    Vector<DType>& scale() {return scale_;}
    const Matrix<DType>& shear() const {return shear_;}
    Matrix<DType>& shear() {return shear_;}

    Matrix<DType> linear() const {return rotation_.apply( mxm::diagonalMatrix(scale_).matmul(shear_));}


private:
    Vector<DType> translation_;
    Matrix<DType> shear_;
    Vector<DType> scale_;
    Rotation<DType, DIM> rotation_;
public:

};


} // namespace mxm



#endif // __TRANSFORM_TRS_H__
