#if !defined(_TEST_ROTATION_H_)
#define _TEST_ROTATION_H_

#include "test_config.h"

#if TEST_AVAILABLE_ALL

#include "mxm/rotation.h"
#include "mxm/string.h"
#include "mxm/lie_special_orthogonal.h"
#include <iostream>

using namespace mxm;

inline void rotationTestCase1()
{
    auto r1 = Rotation<float>::fromAxisAngle(Vec({0,0,1}), 0.4);
    auto r2 = Rotation<float>::fromAxisAngle(Vec({0,0,1}), 0.2);
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
    auto r = Rotation<float>::fromMatrix(mxm::bivectorToRotationMatrix(u1, u2));

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

    auto r(Rotation<float>::fromPlaneAngle(Vec({0,1,0}), Vec({0,0,1}), M_PI));
    if((r.asMatrix() - expect).norm() > 5*eps())
    {
        std::cout << "\n" << mxm::to_string(r.asMatrix());
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void rotationTestCase5()
{
    Matrix<double> rot({3,3},
        {0.999999999999995, -9.999999999999984e-08 ,  0.,
        9.999999999999984e-08 , 0.999999999999995,  0.,
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
        {0.5403023058681398, -0.8414709848078965 ,  0.,
        0.8414709848078965 , 0.5403023058681398,  0.,
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
        auto r(Rotation<float>::fromAxisAngle(Vec({1.,1,1}), 0.1 * i));
        if(fabs(mxm::det(r.asMatrix()) - 1) > 4 * eps())
        {
            std::cout << mxm::det(r.asMatrix()) - 1 << std::endl;
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

    if(1){ // test so3 slerp
        auto mid_expected = rodrigues3D(Vector<float>({1,1,1}), float(M_PI_2 * 0.5));
        auto r1 = rodrigues3D(Vector<float>({1,1,1}), float(M_PI_2));
        auto r0 = Matrix<float>::identity(3);
        auto mid = SO::interp(r0, r1, 0.5f);
        if(norm(mid - mid_expected) > 10 * std::numeric_limits<float>::epsilon())
        {
            std::cout << mxm::to_string(r1) << std::endl;
            std::cout << mxm::to_string(so::exp<3>(1.f * SO::log<3>(r1))) << std::endl;
            std::cout << mxm::to_string(mid) << std::endl;
            std::cout << mxm::to_string(mid_expected) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {
        Quaternion<float> q1({0.407854,0.549993,-0.616159,0.389244});
        Matrix<float> expected_so3({3,3},
            {-0.06232562, -0.99527573, -0.07444288,
            -0.36025683,  0.09199361, -0.92830609,
            0.93076879, -0.0310387 , -0.36428844});

        auto mat_so3 = toSO3(q1);
        if(norm(mat_so3 - expected_so3) > 3 * std::numeric_limits<float>::epsilon())
        {
            std::cout << (norm(mat_so3 - expected_so3)) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    {
        double theta1 = -M_PI + 0.1;
        double theta2 = M_PI - 0.1;

        if(norm(angularDistance(theta1, theta2) + 0.2) > std::numeric_limits<double>::epsilon() * 5)
        {
            std::cout << angularDistance(theta1, theta2) << std::endl;
            std::cout << norm(angularDistance(theta1, theta2) - 0.2) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(norm(angularDistance(theta2, theta1) - 0.2) > std::numeric_limits<double>::epsilon() * 5)
        {
            std::cout << angularDistance(theta2, theta1) << std::endl;
            std::cout << norm(angularDistance(theta2, theta1) - 0.2) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    if(0){
        size_t n = 5;
        Mat rot({n,n}, {
         8.52072892e-01, -1.86443589e-01, -1.52868423e-01, 4.48402930e-01, -1.21559174e-01,
         1.66260826e-01, -3.10093307e-01,  4.89428897e-01, -6.23694169e-02,  7.95467717e-01,
        -2.09350679e-01,  3.08823048e-01,  6.01299790e-01, 6.90088793e-01, -1.51712355e-01,
         4.81917029e-04,  6.33156375e-01, -4.73611309e-01, 2.52124640e-01,  5.57887325e-01,
         4.50001318e-01,  6.10591729e-01,  3.88872076e-01, -5.05228157e-01, -1.34905788e-01});

        auto angles = SO::findAngles(rot);
        std::cout << "angles size: " << angles.size() << std::endl;
        std::cout << mxm::to_string(angles.T()) << std::endl;
    }

    {
        float angle = 0.4;
        auto r1 = Rotation<float>::fromAxisAngle(Vec({0,0,1}), angle);
        auto angles = SO::findAngles(r1.asMatrix());
        if(abs(angles(0) - angle) > std::numeric_limits<float>::epsilon())
        {
            std::cout << "angles size: " << angles.size() << std::endl;
            std::cout << mxm::to_string(angles.T()) << std::endl;
        }

    }

    {
        auto r = Rotation<float, 4>::identity();
    }

#if 0
// special case, single rotation of SO(5)
    {
        size_t n = 5;
        Mat rot({n,n}, {
        0.866734,-0.031498,-0.092792,0.463072,-0.157269,
        0.211789, 0.854830, 0.318394, -0.321836, -0.139493,
        -0.305797, 0.098652, 0.637423, 0.700036, -0.019918,
        -0.052953, 0.401572, -0.399420, 0.305695, 0.763512,
        0.328027, -0.311899, 0.569360, -0.313824, 0.610296});

        std::cout << mxm::to_string(so_n::toAngles(rot)) << std::endl;

        auto skew = so_n::log(rot);
        auto angles = so_n::anglesFromSkewSymmetric(skew);
        std::cout << mxm::to_string(angles) << std::endl;

        auto dianna_mat = so_n::diannaExpMatrix(angles);
        std::cout << "dianna_mat: \n" << mxm::to_string(dianna_mat) << std::endl;
        auto dianna_coeff = so_n::diannaExpCoeff(skew);
        std::cout << "dianna_coeff: \n" << mxm::to_string(dianna_coeff) << std::endl;
        std::cout << "skew: \n" << mxm::to_string(skew) << std::endl;
        auto result =  qr::solve(dianna_mat, dianna_coeff);

        auto b0 = so_n::unsignedWedge(result(Row(0)).T());
        auto b1 = so_n::unsignedWedge(result(Row(1)).T());
        std::cout << "sum: \n" << mxm::to_string(angles(0) * b0 + angles(1) * b1) << std::endl;
        std::cout << "b0: \n" << mxm::to_string(b0) << std::endl;
        std::cout << "b1: \n" << mxm::to_string(b1) << std::endl;

        // std::cout << mxm::to_string(so_n::exp(skew)) << std::endl;
    }
#endif

}
#else
inline void testRotation(){}
#endif
#endif // _TEST_ROTATION_H_
