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
    auto ur5e = ElbowManipulator<double>::buy(eUR5e);
    // auto angles = Vector<double>{-30 * M_PI /180,0,0,0,0,30 * M_PI /180};
    auto angles = Vector<double>{0,0,0,0,0,0};
    auto end_pose = ur5e.forwardKinematics(angles);
    std::vector<Matrix<double>> accum_below;
    std::vector<Matrix<double>> accum_above;
    auto tfs = ur5e.transforms(angles);
    accumulateMatMul(tfs, &accum_below, &accum_above);

    double error(0);
    for(size_t i = 0; i < 6; i++)
    {
        auto local_tf = accum_below.at(i).matmul(tfs.at(i)).matmul(accum_above.at(i));
        if(!isZero(local_tf - end_pose, &error, 50*eps<double>()))
        {
            std::cout << mxm::to_string(local_tf) << std::endl;
            std::cout << "error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

}

void testJacobian()
{
    auto ur5e = ElbowManipulator<double>::buy(eUR5e);
    auto angles = Vector<double>{0,0,0,0,0,0};
    auto end_pose = ur5e.forwardKinematics(angles);

    double inc = 1e-5;
    auto basis = Matrix<double>::identity(6);
    auto jac_num = Matrix<double>::identity(6);
    auto jac_analytical = ur5e.jacob(angles, end_pose);
    for(size_t i = 0; i < 6; i++)
    {
        auto angles_inc = angles + basis(Col(i)) * inc;
        auto end_pose_inc = ur5e.forwardKinematics(angles_inc);
        auto distance = se::vee(SE::log(end_pose_inc.matmul(SE::inv(end_pose))));
        jac_num(Col(i)) = distance / inc;
    }
    double error(0);
    if(!isZero(jac_num - jac_analytical, &error, 1e-7))
    {
        std::cout << mxm::to_string(jac_num) << std::endl;
        std::cout << mxm::to_string(jac_analytical) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testInverseKinematics()
{
    auto ur5e = ElbowManipulator<double>::buy(eUR5e);
    auto angles = Vector<double>{0.1,0.2,0,0,0,0};
    auto guess = Vector<double>{0.11,0.19,0,0,0,0};
    auto desire = ur5e.forwardKinematics(angles);
    auto result = ur5e.inverseKinematics(desire, guess);
    std::cout << mxm::to_string(result) << std::endl;
}

void testJoint()
{
    testUR5e();
    testAccumulativeMatmul();
    testJacobian();
    // testInverseKinematics();
}