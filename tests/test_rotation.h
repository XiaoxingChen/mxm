#if !defined(_TEST_ROTATION_H_)
#define _TEST_ROTATION_H_

#include "mxm/rotation.h"
#include <iostream>

using namespace mxm;

inline void rotationTestCase1()
{
    Rotation r1 = Rotation::fromAxisAngle(Vec({0,0,1}), 0.4);
    Rotation r2 = Rotation::fromAxisAngle(Vec({0,0,1}), 0.2);
    Vec v({1,0,0});
    auto result = (r1 * r2).apply(v);
    Vec expect({0.82533561, 0.56464247, 0.});
    if((result - expect).norm() > eps())
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void rotationTestCase2()
{
    Vec u1({1, 0, 0.1});
    Vec u2({1, 1, 0.1});
    auto r = Rotation::fromMatrix(mxm::bivectorToRotationMatrix(u1, u2));

    Vec v({1,0,1});
    Vec v1(r.apply(v));

    Vec expect({-0.08369046, 1.09452736, 0.89163095});
    if((v1 - expect).norm() > 10*eps())
    {
        // std::cout << mxm::to_string(r.asMatrix()) << std::endl;
        std::cout << mxm::to_string(v1.T()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void rotationTestCase4()
{
    Mat expect(Mat::identity(3));
    expect(1,1) = -1;
    expect(2,2) = -1;

    Rotation r(Rotation::fromPlaneAngle(Vec({0,1,0}), Vec({0,0,1}), M_PI));
    if((r.asMatrix() - expect).norm() > 5*eps())
    {
        std::cout << "\n" << mxm::to_string(r.asMatrix());
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void rotationTestCase5()
{
    Matrix<double> rot({3,3},
        {0.999999999999995, 9.999999999999984e-08 ,  0.,
        -9.999999999999984e-08 , 0.999999999999995,  0.,
        0.        ,  0.        ,  1.});

    Vector<double> axis;
    double angle;
    matrixToAxisAngle3D(rot, axis, angle);
    if((axis - Vector<double>({0,0,1})).norm() > eps() || fabs(angle - 1e-7) > eps())
    {
        std::cout << "axis: " << mxm::to_string(axis.T()) << "angle: " << angle << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void rotationTestCase6()
{
    Matrix<double> rot({3,3},
        {0.5403023058681398, 0.8414709848078965 ,  0.,
        -0.8414709848078965 , 0.5403023058681398,  0.,
        0.        ,  0.        ,  1.});

    Vector<double> axis;
    double angle;
    matrixToAxisAngle3D(rot, axis, angle);
    if((axis - Vector<double>({0,0,1})).norm() > eps() || fabs(angle - 1) > eps())
    {
        std::cout << "axis: " << mxm::to_string(axis.T()) << "angle: " << angle << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void rotationTestDeterminant()
{
    for(size_t i = 0; i < 10; i++)
    {
        Rotation r(Rotation::fromAxisAngle(Vec({1.,1,1}), 0.1 * i));
        if(fabs(r.asMatrix().det() - 1) > 4 * eps())
        {
            std::cout << r.asMatrix().det() - 1 << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }
}

inline void testRotation()
{
    rotationTestCase1();
    rotationTestCase2();
    rotationTestCase4();
    rotationTestCase5();
    rotationTestCase6();
    rotationTestDeterminant();
}

#endif // _TEST_ROTATION_H_
