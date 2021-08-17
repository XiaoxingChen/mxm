#include "mxm/kinematics_manipulator.h"
using namespace mxm;

void testUR5e()
{
    auto ur5e = ElbowManipulator<float>::buy(eUR5e);
    auto joint_degrees = Vector<float>{30, 30,30,30,30,30};
    auto rot_dir = Vector<float>{-1, -1,-1,-1,-1,-1};
    auto joint_angles = rot_dir * joint_degrees * M_PI / 180;
    auto end_pose = ur5e.forwardKinematics(joint_angles);
    auto expected_pose = Matrix<float>({4,4},
        {-0.21651, -0.87500, 0.43301, -0.29246,
        -0.62500, -0.21651, -0.75000, -0.42237,
        0.75000, -0.43301, -0.50000, -0.43946,
        0.00000, 0.00000, 0.00000, 1.00000});


    float error(0);
    if(!isZero(end_pose- expected_pose, &error, 50*eps<float>()))
    {
        std::cout << "tfs: " << mxm::to_string(ur5e.transforms(joint_angles)) << std::endl;
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
    auto angles = Vector<double>{0.,0.,M_PI_2,0.,0.5,0.6};
    double rot_dir = -1.;
    auto end_pose = ur5e.forwardKinematics(angles * rot_dir);

    Matrix<double> dhParams({6, dh::PARAM_NUM}, {
            0, 0, 0.1625, M_PI_2,
            0, -0.425, 0, 0,
            0, -0.3922, 0, 0,
            0, 0, 0.1333, M_PI_2,
            0, 0, 0.0997, -M_PI_2,
            0, 0, 0.0996, 0}, ROW);

    auto angle0 = ur5::ikAngle0(dhParams, end_pose);
    std::cout << "angle0: " << mxm::to_string(angle0) << std::endl;
    dhParams(0, dh::IDX_THETA) = angle0(0);
    auto angle4 = ur5::ikAngle4(dhParams, end_pose);
    std::cout << "angle4: " << mxm::to_string(angle4) << std::endl;

    dhParams(4, dh::IDX_THETA) = angle4(0);

    auto angle5 = ur5::ikAngle5(dhParams, end_pose);
    std::cout << "angle5: " << (angle5) << std::endl;
    dhParams(5, dh::IDX_THETA) = angle5;

    auto angle2 = ur5::ikAngle2(dhParams, end_pose);
    std::cout << "angle2: " << mxm::to_string(angle2) << std::endl;

    dhParams(2, dh::IDX_THETA) = angle2(0);

    auto angle1 = ur5::ikAngle1(dhParams, end_pose);
    std::cout << "angle1: " << (angle1) << std::endl;
    dhParams(1, dh::IDX_THETA) = angle1;


}

void testJoint()
{
    testUR5e();
    testAccumulativeMatmul();
    testJacobian();
    // testInverseKinematics();
}