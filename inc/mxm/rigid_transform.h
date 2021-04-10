#if !defined(_RIGID_TRANSFORM_H_)
#define _RIGID_TRANSFORM_H_

#include "linalg.h"
#include "full_dimensional_rotation.h"

namespace mxm
{
template<typename DType>
class RigidTransform
{
public:
    using Rotation = FullDimensionalRotation<DType>;
    using Translation = Vector<DType>;
    RigidTransform(const Translation& p, const Rotation& r): translation_(p), rotation_(r)
    {
        if(p.size() != r.dim())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    RigidTransform(const Rotation& r): translation_(Vector<DType>::zeros(r.dim())), rotation_(r)
    {}

    RigidTransform(const Translation& p): translation_(p), rotation_(Rotation::Identity(p.size()))
    {}

    using ThisType = RigidTransform;

    ThisType operator*(const ThisType& rhs) const { return ThisType(apply(rhs.translation_) , rotation_ * rhs.rotation_); }

    Matrix<DType> apply(const Matrix<DType>& vectors) const
    {
        auto ret = rotation_.apply(vectors);
        for(size_t i = 0; i < vectors.shape(1); i++)
            ret(Col(i)) += translation_;
        return ret;
    }

    Matrix<DType> asMatrix() const
    {
        Matrix<DType> ret = Matrix<DType>::Identity(dim() + 1);
        ret.setBlock(0,0, rotation_.asMatrix());
        ret.setBlock(0, dim(), translation_);
        return ret;
    }

    ThisType inv() const { return ThisType(rotation_.inv().apply(-translation_), rotation_.inv()); }

    size_t dim() const { return translation_.size(); }

    static ThisType Identity(size_t dim) { return ThisType(Translation::zeros(dim), Rotation::Identity(dim)); }

    const Rotation& rotation() const {return rotation_;}
    Rotation& rotation() {return rotation_;}
    const Translation& translation() const {return translation_;}
    Translation& translation() {return translation_;}

private:
    Translation translation_;
    Rotation rotation_;
};

using RigidTrans = RigidTransform<FloatType>;
} // namespace mxm


#endif // _RIGID_TRANSFORM_H_
