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
        // std::cout << r.asMatrix().str() << std::endl;
        std::cout << v1.T().str() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void rotationTestCase4()
{
    Mat expect(Mat::Identity(3));
    expect(1,1) = -1;
    expect(2,2) = -1;

    Rotation r(Rotation::fromPlaneAngle(Vec({0,1,0}), Vec({0,0,1}), M_PI));
    if((r.asMatrix() - expect).norm() > 5*eps())
    {
        std::cout << "\n" << r.asMatrix().str();
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void rotationTestDeterminant()
{
    for(size_t i = 0; i < 10; i++)
    {
        Rotation r(Rotation::fromAxisAngle(UnitVec({1.,1,1}), 0.1 * i));
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
    rotationTestDeterminant();
}

#endif // _TEST_ROTATION_H_
