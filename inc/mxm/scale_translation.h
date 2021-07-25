#if !defined(_SCALE_TRANSLATION_)
#define _SCALE_TRANSLATION_

#include "linalg.h"

namespace mxm
{

template <typename DType>
class ScaleTranslation
{
private:
    Vector<DType> s_;
    Vector<DType> t_;
public:
    ScaleTranslation(const Vector<DType>& scale, const Vector<DType>& translation);

    using ThisType = ScaleTranslation<DType>;

    const Vector<DType>& scale() const { return s_; }
    Vector<DType>& scale() { return s_; }

    const Vector<DType>& translation() const { return t_; }
    Vector<DType>& translation() { return t_; }

    Matrix<DType> apply(const Matrix<DType>& vecs) const;

    ThisType inv() const;
    Matrix<DType> asMatrix() const;
    static Matrix<DType> inv(const Matrix<DType>& mat);
    static ThisType fromMatrix(const Matrix<DType>& mat);
    static ThisType fromTranslationBeforeScale(const Vector<DType>& translation, const Vector<DType>& scale);
};

template <typename DType>
ScaleTranslation<DType>::ScaleTranslation(const Vector<DType>& scale, const Vector<DType>& translation)
:s_(scale), t_(translation)
{
    if(1 == scale_.size())
    {
        s_ = Vector<DType>::ones(translation.size());
        s_ *= scale;
    }
    assert(s_.size() == t_.size());
}

template <typename DType>
Matrix<DType> ScaleTranslation<DType>::apply(const Matrix<DType>& vecs) const
{
    assert(vecs.shape(0) <= s_.size())
    if(vecs.shape(0) == s_.size())
    {
        return vecs * s_ + t_;
    }
    return vecs(Block({0, s_.size()}, {})) * s_ + t_;
}

template <typename DType>
ScaleTranslation<DType> ScaleTranslation<DType>::inv() const
{
    Vector<DType> inv_scale = DType(1.) / scale();
    Vector<DType> inv_trans = translation() * inv_scale * DType(-1.);
    return ThisType(inv_scale, inv_trans);
}

template <typename DType>
Matrix<DType> ScaleTranslation<DType>::inv(const Matrix<DType>& mat)
{
    assert(mat.square());
    size_t n = mat.shape(0);
    Matrix<DType> ret = Matrix<DType>::identity(n);
    for(size_t i = 0; i < n; i++)
    {
        ret(i,i) = DType(1.) / mat(i,i);
        ret(i,n-1) = DType(-1.) * mat(i, n-1) / mat(i,i);
    }
    return ret;
}

template <typename DType>
Matrix<DType> ScaleTranslation<DType>::asMatrix() const
{
    size_t n = s_.shape(0);
    Matrix<DType> ret = Matrix<DType>::identity(n);
    for(size_t i = 0; i < n; i++)
    {
        ret(i,i) = s_(i);
        ret(i,n-1) = t_(i);
    }
    return ret;
}

template <typename DType>
ScaleTranslation<DType> ScaleTranslation<DType>::fromMatrix(const Matrix<DType>& mat)
{
    assert(false && "not implemented yet");
    Vector<DType> scale;
    Vector<DType> translation;
    return ThisType(scale, translation);
}

template <typename DType> ScaleTranslation<DType>
ScaleTranslation<DType>::fromTranslationBeforeScale(
    const Vector<DType>& translation,
    const Vector<DType>& scale)
{
    return ThisType(scale, Vector<DType>(translation / scale));
}


} // namespace mxm

#endif // _SCALE_TRANSLATION_
