#include "mxm/kinematics_manipulator.h"
using namespace mxm;

void testUR5e()
{
    auto ur5e = ElbowManipulator<float>::buy(eUR5e);
    auto end_pose = ur5e.forwardKinematics({-30 * M_PI /180,0,0,0,0,30 * M_PI /180});
    auto expected_pose = Matrix<float>({4,4},
        {0.75,-0.43301,0.5,-0.59127,
        0.43301, -0.25000, -0.86603, -0.61030,
        0.50000, 0.86603, 0.00000, 0.06280,
        0.00000, 0.00000, 0.00000, 1.00000});


    float error(0);
    if(!isZero(end_pose- expected_pose, &error, 50*eps<float>()))
    {
        std::cout << mxm::to_string(end_pose) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testAccumulativeMatmul()
{
    auto ur5e = ElbowManipulator<float>::buy(eUR5e);
    auto angles = Vector<float>{-30 * M_PI /180,0,0,0,0,30 * M_PI /180};
    auto end_pose = ur5e.forwardKinematics(angles);
    std::vector<Matrix<float>> accum_below;
    std::vector<Matrix<float>> accum_above;
    auto tfs = ur5e.transforms(angles);
    accumulateMatMul(tfs, &accum_below, &accum_above);

    float error(0);
    for(size_t i = 0; i < 6; i++)
    {
        auto local_tf = accum_below.at(i).matmul(tfs.at(i)).matmul(accum_above.at(i));
        if(!isZero(local_tf - end_pose, &error, 50*eps<float>()))
        {
            std::cout << mxm::to_string(local_tf) << std::endl;
            std::cout << "error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    // std::cout << mxm::to_string(ur5e.jacob(angles, end_pose)) << std::endl;

    // std::cout << mxm::to_string(accum_above) << std::endl;

}

void testJoint()
{
    testUR5e();
    testAccumulativeMatmul();
}