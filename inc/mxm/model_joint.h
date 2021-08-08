#if !defined(__MODEL_JOINT_H__)
#define __MODEL_JOINT_H__

#include "rigid_transform.h"
#include "rotation.h"
#include <bitset>

namespace mxm
{
template<typename DType>
class Joint
{
public:

    using ThisType = Joint<DType>;
    virtual void set(FloatType val) = 0;
    Joint(const RigidTransform<DType>& zero_tf=RigidTransform<DType>::identity(3), bool install_dir=true):
        zero_tf_(zero_tf), state_tf_(RigidTransform<DType>::identity(3)), inverse_install_(!install_dir){}

    RigidTransform<DType> tf() const
    {
        auto ret = zero_tf_ * state_tf_;
        if(inverse_install_) return ret.inv();
        return ret;
    }

    void setAngle(DType val) { state_tf_.rotation() = Rotation<DType>::fromAxisAngle(Vector<DType>{0,0,1}, val); }
    RigidTransform<DType> zeroPosition() const
    {
        if(inverse_install_) return zero_tf_.inv();
        return zero_tf_;
    }

    bool inverseInstall() const { return inverse_install_; }

protected:
    RigidTransform<DType> zero_tf_;
    RigidTransform<DType> state_tf_;
    bool inverse_install_;
};


} // namespace mxm



#endif // __MODEL_JOINT_H__
