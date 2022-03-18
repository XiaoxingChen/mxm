#if !defined(_RIGID_TRANSFORM_H_)
#define _RIGID_TRANSFORM_H_

#include "linalg.h"
#include "full_dimensional_rotation.h"

namespace mxm
{
template<typename DType, size_t DIM=3>
class RigidTransform
{
public:

    using Translation = Vector<DType>;
    RigidTransform(const Translation& p, const Rotation<DType, DIM>& r): translation_(p), rotation_(r)
    {
        if(p.size() != r.dim())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    RigidTransform(const Rotation<DType, DIM>& r): translation_(Vector<DType>::zeros(r.dim())), rotation_(r)
    {}

    RigidTransform(const Translation& p): translation_(p), rotation_(Rotation<DType, DIM>::identity())
    {}

    template<typename RhsDType=DType>
    RigidTransform(const RigidTransform<RhsDType, DIM>& rhs)
    {
        translation_ = rhs.translation();
        rotation_ = rhs.rotation();
    }

    template<typename RhsDType=DType>
    void operator = (const RigidTransform<RhsDType, DIM>& rhs)
    {
        translation_ = rhs.translation();
        rotation_ = rhs.rotation();
    }


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
        Matrix<DType> ret = Matrix<DType>::identity(dim() + 1);
        ret.setBlock(0,0, rotation_.asMatrix());
        ret.setBlock(0, DIM, translation_);
        return ret;
    }

    ThisType inv() const { return ThisType(rotation_.inv().apply(-translation_), rotation_.inv()); }

    constexpr size_t dim() const { return DIM; }

    static ThisType identity() { return ThisType(Translation::zeros(DIM), Rotation<DType, DIM>::identity()); }
    static ThisType fromMatrix(const Matrix<DType>& se)
    {
        auto tra = Vector<DType>(se(Block({0, DIM}, {DIM, DIM + 1})));
        auto rot = Rotation<DType, DIM>::fromMatrix(se(Block({0, DIM}, {0, DIM})));
        return ThisType(tra, rot);
    }

    const Rotation<DType, DIM>& rotation() const {return rotation_;}
    Rotation<DType, DIM>& rotation() {return rotation_;}
    const Translation& translation() const {return translation_;}
    Translation& translation() {return translation_;}

private:
    Translation translation_;
    Rotation<DType, DIM> rotation_;
};

using RigidTrans = RigidTransform<FloatType>;
} // namespace mxm


#endif // _RIGID_TRANSFORM_H_
