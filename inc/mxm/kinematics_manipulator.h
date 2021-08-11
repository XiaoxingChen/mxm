#if !defined(_KINEMATICS_MANIPULATOR_H_)
#define _KINEMATICS_MANIPULATOR_H_

#include "linalg.h"
#include "rigid_transform.h"
#include "lie_special_euclidean.h"
#include <bitset>


namespace mxm
{

enum ManipulatorModels
{
    eUR3e = 0,
    eUR5e,
    eUR10e
};

template <typename DType>
void accumulateMatMul(
    const std::vector<Matrix<DType>>& tfs,
    std::vector<Matrix<DType>>* accum_below,
    std::vector<Matrix<DType>>* accum_above)
{
    size_t N = tfs.front().shape(0);
    accum_below->resize(tfs.size());
    accum_above->resize(tfs.size());

    accum_above->back() = Matrix<DType>::identity(N);
    accum_below->front() = Matrix<DType>::identity(N);
    for(size_t i = 1; i < tfs.size(); i++)
    {
        accum_below->at(i) = accum_below->at(i-1).matmul(tfs.at(i-1));
        size_t r_i = tfs.size() - 1 - i;
        accum_above->at(r_i) = tfs.at(r_i+1).matmul(accum_above->at(r_i+1));
    }
}

// reference:
// docs/derivative_of_manipulator_se3.html
template <typename DType> Matrix<DType>
jointJacobian(
    const Matrix<DType>& accum_below,
    const Matrix<DType>& accum_above,
    const Matrix<DType>& zero_pose,
    const Matrix<DType>& axis_rot,
    const Matrix<DType>& inv_desire,
    bool axis_out)
{
    if(axis_out)
    {
        auto term1 = accum_below.matmul(zero_pose);
        auto term2 = accum_above.matmul(inv_desire);
        auto term_log = term1.matmul(axis_rot).matmul(term2);
        return se::jacobInv<3>(SE::log<3>(term_log)).matmul(SE::adj<3>(term1));
    }
    auto inv_rot = SE::inv<3>(axis_rot);
    auto term1 = accum_below;
    auto term2 = SE::inv<3>(zero_pose).matmul(accum_above).matmul(inv_desire);
    auto term_log = term1.matmul(inv_rot).matmul(term2);
    return -se::jacobInv<3>(SE::log<3>(term_log)).matmul(SE::adj<3>(term1.matmul(inv_rot)));
}

// https://www.universal-robots.com/articles/ur/application-installation/dh-parameters-for-calculations-of-kinematics-and-dynamics/
template<typename DType=double>
class ElbowManipulator
{
public:
    static const size_t JOINT_NUM = 6;
    using ThisType = ElbowManipulator<DType>;
    ElbowManipulator(std::vector<Matrix<DType>>&& zero_pose, const std::bitset<6>& axis_out)
        :axis_out_(axis_out)
    {
        zero_pose_.reserve(JOINT_NUM);
        for(size_t i = 0; i < JOINT_NUM; i++)
        {
            zero_pose_.push_back(RigidTransform<DType>::fromMatrix(zero_pose.at(i)));
        }
    }

    std::vector<Matrix<DType>> transforms(const Vector<DType>& angles) const
    {
        std::vector<Matrix<DType>> ret(JOINT_NUM);
        Vector<DType> z_axis{0,0,1};
        for(size_t i = 0; i < JOINT_NUM; i++)
        {
            auto angle_pose = RigidTransform<DType>(Rotation<DType>::fromAxisAngle(z_axis, angles(i)));
            RigidTransform<DType> local_pose = zero_pose_.at(i) * angle_pose;
            if(axis_out_[i])
                ret.at(i) = local_pose.asMatrix();
            else
                ret.at(i) = local_pose.inv().asMatrix();
        }
        return ret;
    }

    Matrix<DType> forwardKinematics(const Vector<DType>& angles) const
    {
        Matrix<DType> ret = Matrix<DType>::identity(4);
        auto tfs = transforms(angles);
        for(size_t i = 0; i < JOINT_NUM; i++)
        {
            ret = ret.matmul(tfs.at(i));
        }
        return ret;
    }

    Matrix<DType> jacob(const Vector<DType>& angles, const Matrix<DType>& desire) const
    {
        Matrix<DType> jac({6, JOINT_NUM});
        std::vector<Matrix<DType>> accum_below;
        std::vector<Matrix<DType>> accum_above;
        auto inv_desire = SE::inv<3>(desire);
        auto tfs = transforms(angles);
        accumulateMatMul(tfs, &accum_below, &accum_above);
        Vector<DType> z_axis{0,0,1};
        for (size_t i = 0; i < JOINT_NUM; i++)
        {
            auto zero_pose = zero_pose_.at(i).asMatrix();
            // auto angle_pose = SE::inv(zero_pose).matmul(tfs.at(i));
            auto angle_pose = RigidTransform<DType>(Rotation<DType>::fromAxisAngle(z_axis, angles(i))).asMatrix();
            auto joint_jac = jointJacobian(accum_below.at(i), accum_above.at(i),zero_pose, angle_pose, inv_desire, axis_out_[i]);
            jac(Col(i)) = joint_jac(Col(end() - 1));
        }

        return jac;
    }

    // std::array<

    // Vector<DType> inverseKinematics(const Matrix<DType>& se3) const
    // {

    // }
    static ThisType buy(const ManipulatorModels&);

private:
    std::vector<RigidTransform<DType>> zero_pose_;
    std::bitset<JOINT_NUM> axis_out_;
};

// Reference: https://en.wikipedia.org/wiki/Denavit%E2%80%93Hartenberg_parameters#Denavit%E2%80%93Hartenberg_matrix
template<typename DType>
std::vector<Matrix<DType>> denavitHartenbergMatrix(const Matrix<DType>& params)
{
    std::vector<Matrix<DType>> ret(params.shape(1));
    for(size_t j = 0; j < params.shape(1); j++)
    {
        const auto & theta = params(0, j);
        const auto & r     = params(1, j);
        const auto & d     = params(2, j);
        const auto & alpha = params(3, j);
        ret.at(j) = Matrix<DType>({4,4},{
            cos(theta), -sin(theta) * cos(alpha), sin(theta) * sin(alpha), r * cos(theta),
            sin(theta), cos(theta) * cos(alpha), -cos(theta) * sin(alpha), r * sin(theta),
            0         , sin(alpha)              , cos(alpha)             , d,
            0, 0, 0, 1});
    }
    return ret;
}

template<typename DType>
ElbowManipulator<DType> ElbowManipulator<DType>::buy(const ManipulatorModels& model)
{

    std::vector<Matrix<DType>> zero_poses(JOINT_NUM);
    if(eUR5e == model)
    {
        static const Matrix<DType> dhParams({4,JOINT_NUM}, {
            0, 0, 0.1625, M_PI_2,
            0, -0.425, 0, 0,
            0, -0.3922, 0, 0,
            0, 0, 0.1333, M_PI_2,
            0, 0, 0.0997, -M_PI_2,
            0, 0, 0.0996, 0}, COL);

        auto zero_poses = denavitHartenbergMatrix(dhParams);
        std::bitset<JOINT_NUM> inversion(0x3);
        for(size_t i = 0; i < 6; i++)
        {
            if(inversion[i])
            {
                zero_poses.at(i) = SE::inv(zero_poses.at(i));
            }
        }
        return ThisType(std::move(zero_poses), ~inversion);
    }
    return ThisType(std::move(zero_poses), std::bitset<JOINT_NUM>(0));
}

} // namespace mxm

#endif // _KINEMATICS_MANIPULATOR_H_
