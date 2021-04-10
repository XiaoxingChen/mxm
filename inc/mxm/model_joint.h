#if !defined(__MODEL_JOINT_H__)
#define __MODEL_JOINT_H__

#include "rigid_transform.h"

namespace mxm
{
class Joint
{
public:
    virtual void set(FloatType val) = 0;
    Joint(const RigidTrans& zero_tf, bool rotor_output=true):
        zero_tf_(zero_tf), state_tf_(RigidTrans::Identity(3)), rotor_output_(rotor_output){}

    RigidTrans nominalStateTf() const
    {
        if(rotor_output_) return state_tf_;
        return state_tf_.inv();
    }

    RigidTrans output() const
    {
        return zero_tf_ * nominalStateTf();
    }

protected:
    RigidTrans zero_tf_;
    RigidTrans state_tf_;
    bool rotor_output_;
};

// rotate around z-axis
class AngularJoint: public Joint
{
public:
    virtual void set(FloatType val) override
    {
        state_tf_.rotation() = RigidTrans::Rotation::fromAxisAngle(Vec({0,0,1}), val);
    }
private:

};
} // namespace mxm



#endif // __MODEL_JOINT_H__
