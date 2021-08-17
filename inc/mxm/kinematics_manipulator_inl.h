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
    virtual Vector<DType> residualOnPerturb(const Vector<DType>& increment) const override
    {
        auto perturbed_angle = joint_angles_ + increment;
        auto end_pose = manipulator_.forwardKinematics(perturbed_angle);
        return se::vee(SE::log<3>(end_pose.matmul(inv_desire_)));
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

namespace ur5
{
// Be careful that the joint indices start form 0, which is different from the reference paper.
// reference:
// [1] http://rasmusan.blog.aau.dk/files/ur5_kinematics.pdf
template<typename DType>
Vector<DType> ikAngle0(const Matrix<DType>& dh_param, const Matrix<DType>& desire)
{
    Vector<DType> solutions = Vector<DType>::zeros(2);
    // formula (5)
    Vector<DType> p4 = desire.matmul(Vector<DType>{0, 0, -dh_param(5, dh::IDX_D), 1});

    // formula (9)
    solutions += (M_PI_2 + atan2(p4(1), p4(0)));
    auto phi = acos(dh_param(3, dh::IDX_D) / sqrt(p4(0) * p4(0) + p4(1) * p4(1)));
    solutions(0) += phi;
    solutions(1) -= phi;
    return solutions;
}

template Vector<float> ikAngle0(const Matrix<float>& dh_param, const Matrix<float>& desire);
template Vector<double> ikAngle0(const Matrix<double>& dh_param, const Matrix<double>& desire);

template<typename DType>
Vector<DType> ikAngle4(const Matrix<DType>& dh_param, const Matrix<DType>& desire)
{
    Vector<DType> solutions = Vector<DType>::zeros(2);
    Vector<DType> p5 = desire(Col(3));

    const DType& angle0 = dh_param(0, dh::IDX_THETA);
    // formula (12)
    DType propotion = (p5(0) * sin(angle0) - p5(1) * cos(angle0) - dh_param(3, dh::IDX_D)) / dh_param(5, dh::IDX_D);
    solutions(0) = acos(propotion);
    solutions(1) = -acos(propotion);
    return solutions;

}

template Vector<float> ikAngle4(const Matrix<float>& dh_param, const Matrix<float>& desire);
template Vector<double> ikAngle4(const Matrix<double>& dh_param, const Matrix<double>& desire);

template<typename DType>
DType ikAngle5(const Matrix<DType>& dh_param, const Matrix<DType>& desire)
{
    const DType& angle0 = dh_param(0, dh::IDX_THETA);
    const DType& angle4 = dh_param(4, dh::IDX_THETA);
    // formula (16)
    auto desire_inv = SE::inv(desire);
    auto term_y = (-desire_inv(1,0) * sin(angle0) + desire_inv(1,1)*cos(angle0)) / sin(angle4);
    auto term_x = (desire_inv(0,0) * sin(angle0) - desire_inv(0,1)*cos(angle0)) / sin(angle4);
    return atan2(term_y, term_x);
}

template float ikAngle5(const Matrix<float>& dh_param, const Matrix<float>& desire);
template double ikAngle5(const Matrix<double>& dh_param, const Matrix<double>& desire);

template<typename DType>
Vector<DType> ikAngle2(const Matrix<DType>& dh_param, const Matrix<DType>& desire)
{
    auto tfs = denavitHartenbergMatrix(dh_param);
    auto pose_joint3 = desire.matmul(SE::inv(tfs.at(5))).matmul(SE::inv(tfs.at(4)));
    auto p03 = SE::inv(pose_joint3).matmul(tfs.at(0));

    DType p03xz_2 = p03(0, 3) * p03(0, 3) + p03(2, 3) * p03(2, 3);

    Vector<DType> solutions = Vector<DType>::zeros(2);
    const auto& a1 = dh_param(1, dh::IDX_R);
    const auto& a2 = dh_param(2, dh::IDX_R);
    // formula (19)
    auto proportion = (p03xz_2 - a1*a1 - a2*a2) / (2*a1*a2);

    solutions(0) = acos(proportion);
    solutions(1) = -acos(proportion);
    return solutions;
}

template Vector<float> ikAngle2(const Matrix<float>& dh_param, const Matrix<float>& desire);
template Vector<double> ikAngle2(const Matrix<double>& dh_param, const Matrix<double>& desire);

template<typename DType>
DType ikAngle1(const Matrix<DType>& dh_param, const Matrix<DType>& desire)
{
    auto tfs = denavitHartenbergMatrix(dh_param);
    auto pose_joint3 = desire.matmul(SE::inv(tfs.at(5))).matmul(SE::inv(tfs.at(4)));
    auto p03 = SE::inv(pose_joint3).matmul(tfs.at(0));
    // p03 = SE::inv(p03);
    std::cout << "p0:\n" << mxm::to_string(tfs.at(0)) << std::endl;
    std::cout << "p03:\n" << mxm::to_string(p03) << std::endl;

    DType p03xz = sqrt(p03(0, 3) * p03(0, 3) + p03(2, 3) * p03(2, 3));

    const auto& angle2 = dh_param(2, dh::IDX_THETA);
    const auto& a2 = dh_param(2, dh::IDX_R);
    // formula (22)
    auto phi1 = atan2(-p03(2, 3), -p03(0, 3));
    auto phi2 = asin(-a2 * sin(angle2)/p03xz);
    std::cout << "phi1, phi2: " << phi1 << "," << phi2 << std::endl;
    return phi1 - phi2;
}

template float ikAngle1(const Matrix<float>& dh_param, const Matrix<float>& desire);
template double ikAngle1(const Matrix<double>& dh_param, const Matrix<double>& desire);

} // namespace ur5


} // namespace mxm


#endif // _KINEMATICS_MANIPULATOR_INL_H_
