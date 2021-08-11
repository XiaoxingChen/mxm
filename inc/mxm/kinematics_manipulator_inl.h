#if !defined(_KINEMATICS_MANIPULATOR_INL_H_)
#define _KINEMATICS_MANIPULATOR_INL_H_

#if defined(MXM_COMPILED_LIB)
#include "kinematics_manipulator.h"
#endif // MXM_COMPILED_LIB

#include "optimize.h"

namespace mxm
{

template <typename DType>
class ElbowManipulatorInverseKinematics: public opt::NoneLinearProblem<DType>
{
public:
    ElbowManipulatorInverseKinematics(const ElbowManipulator<DType>& manipulator, const Matrix<DType>& desire_end_pose)
    :manipulator_(manipulator), desire_(desire_end_pose), inv_desire_(SE::inv(desire_))
    {
    }
    void initGuess(const Vector<DType>& joint_angles)
    {
        joint_angles_ = joint_angles;
        update(Vector<DType>::zeros(ElbowManipulator<DType>::JOINT_NUM));
    }
    virtual void update(const Vector<DType>& increment) override
    {
        joint_angles_ += increment;
        auto end_pose = manipulator_.forwardKinematics(joint_angles_);
        this->residual_ = se::vee(SE::log<3>(end_pose.matmul(inv_desire_)));
        this->jacobian_ = manipulator_.jacob(joint_angles_, desire_);
    }
    const Vector<DType>& angles() const { return joint_angles_; }
private:
    Vector<DType> joint_angles_;
    ElbowManipulator<DType> manipulator_;
    Matrix<DType> desire_;
    Matrix<DType> inv_desire_;
};

template <typename DType> Matrix<DType>
ElbowManipulator<DType>::inverseKinematics(const Matrix<DType>& desire, const Vector<DType>& guess) const
{
    ElbowManipulatorInverseKinematics<DType> problem(*this, desire);
    problem.initGuess(guess);
    problem.solve(5, 0, "lm");
    return problem.angles();
}

template class ElbowManipulator<double>;
template class ElbowManipulator<float>;

} // namespace mxm


#endif // _KINEMATICS_MANIPULATOR_INL_H_
