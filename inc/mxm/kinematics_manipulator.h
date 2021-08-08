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


    Matrix<DType> forwardKinematics(const Vector<DType>& angles) const
    {
        RigidTransform<DType> ret = RigidTransform<DType>::identity();
        Vector<DType> z_axis{0,0,1};
        for(size_t i = 0; i < JOINT_NUM; i++)
        {
            auto angle_pose = RigidTransform<DType>(Rotation<DType>::fromAxisAngle(z_axis, angles(i)));
            RigidTransform<DType> local_pose = zero_pose_.at(i) * angle_pose;
            if(axis_out_[i])
            {
                ret = ret * local_pose;
            }else
            {
                ret = ret * local_pose.inv();
            }
        }
        return ret.asMatrix();
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

template<typename DType>
ElbowManipulator<DType> ElbowManipulator<DType>::buy(const ManipulatorModels& model)
{

    std::vector<Matrix<DType>> zero_poses(JOINT_NUM);
    if(eUR5e == model)
    {
        std::bitset<JOINT_NUM> inversion(0x3);
        zero_poses.at(0) = Matrix<DType>({4,4},
            {1,0,0,0,
            0,0,-1,0,
            0,1,0,0.1625,
            0,0,0,1});
        zero_poses.at(1) = Matrix<DType>({4,4},
            {1,0,0,-0.425,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1});

        zero_poses.at(2) = Matrix<DType>({4,4},
            {1,0,0,-0.3922,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1});

        zero_poses.at(3) = Matrix<DType>({4,4},
            {1,0,0,0,
            0,0,-1,0,
            0,1,0,0.1333,
            0,0,0,1});
        zero_poses.at(4) = Matrix<DType>({4,4},
            {1,0,0,0,
            0,0,1,0,
            0,-1,0,0.0997,
            0,0,0,1});
        zero_poses.at(5) = Matrix<DType>({4,4},
            {1,0,0,0,
            0,1,0,0,
            0,0,1,0.0996,
            0,0,0,1});
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
