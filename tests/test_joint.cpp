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

void testJoint()
{
    testUR5e();
}