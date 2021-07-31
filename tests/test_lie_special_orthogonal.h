#if !defined(_TEST_LIE_SPECIAL_ORTHOGONAL_H_)
#define _TEST_LIE_SPECIAL_ORTHOGONAL_H_

#include "mxm/lie_special_orthogonal.h"
#include "mxm/lie_special_euclidean.h"

using namespace mxm;

inline Matrix<float> testDataSO3()
{
    // theta = 0.5 * pi
    // axis = [1,1,1]
    return Matrix<float>({3,3},{
    0.333333, -0.244017, 0.910684,
    0.910684, 0.333333, -0.244017,
    -0.244017, 0.910684, 0.333333});
}

inline Matrix<float> testDataSO2()
{
    float c = cos(M_PI * 0.25);
    float s = sin(M_PI * 0.25);
    return Matrix<float>({2,2}, {c, -s, s, c}, ROW);
}

inline Matrix<float> testDataSO5()
{
    size_t n = 5;
    const Matrix<float> rot({n,n}, {
            0.8606430292, -0.0866019130, -0.1694403887, 0.4272221625, -0.2014071941,
            0.2661026716, 0.8530187011, 0.2963767946, -0.3091226518, -0.1347314417,
            -0.2259227931, 0.1157709658, 0.6333417892, 0.7310451865, 0.0026819520,
            -0.0206990279, 0.3864032328, -0.4326551259, 0.3044715226, 0.7552289963,
            0.3701532483, -0.3196074963, 0.5432665348, -0.3078871667, 0.6090193987});

    return rot;
}

inline void testLieSpecialOrthogonal()
{
    {
        auto mat = testDataSO2();
        auto theta = SO::findAngle<2>(mat);
        if(norm(theta - 0.25 * M_PI) > std::numeric_limits<float>::epsilon())
        {
            std::cout << norm(theta - 0.25) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {
        auto rot = testDataSO5();
        auto inv_inv = so::exp(SO::log(rot));
        if(norm(inv_inv - rot) > 15 * std::numeric_limits<float>::epsilon())
        {
            std::cout << mxm::to_string(inv_inv) << std::endl;
            std::cout << "error: " << (norm(inv_inv - rot)) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {
        auto rot = testDataSO3();

        auto inv_inv = so::exp<3>(SO::log<3>(rot));

        if(norm(inv_inv - rot) > 5 * std::numeric_limits<float>::epsilon())
        {
            std::cout << SO::findAngle<3>(rot) << std::endl;
            std::cout << mxm::to_string(SO::log<3>(rot)) << std::endl;
            std::cout << mxm::to_string(inv_inv) << std::endl;
            std::cout << "error: " << (norm(inv_inv - rot)) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    { // test jacob left inv
        using DType = double;
        auto skew1 = so::wedge(Vector<double>{1,2,3} * eps<DType>());
        auto skew2 = so::wedge(Vector<double>{2,1,1});

        auto direct_log = SO::log<3>( so::exp<3>(skew1).matmul(so::exp<3>(skew2)));
        auto jac_inv_log = so::jacobInv(skew2).matmul(skew1) + skew2;



        DType error(0);
        if(!isZero(direct_log - jac_inv_log, &error, 10*eps<DType>()))
        {
            std::cout << mxm::to_string(direct_log) << std::endl;
            std::cout << mxm::to_string(jac_inv_log) << std::endl;
            std::cout << "error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        auto expect_identity = so::jacobInv(skew2).matmul(so::jacob(skew2));

        if(!isIdentity(expect_identity, &error))
        {
            std::cout << "error: " << error << std::endl;
            std::cout << mxm::to_string(expect_identity) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

inline void testLieSpecialEuclidean()
{
    {
        using DType = double;
        auto rot = so::exp<3>(so::wedge(Vector<double>{2,1,1}));
        auto pose = Matrix<DType>::identity(4);
        pose(se::rotBlk<3>()) = rot;
        pose(se::traBlk<3>()) = Vector<DType>{1,2,3};

        auto pose_sym = se::exp<3>(SE::log<3>(pose));

        DType error(0);
        if(!isZero(pose - pose_sym, &error, 10*eps<DType>()))
        {
            std::cout << mxm::to_string(pose) << std::endl;
            std::cout << mxm::to_string(pose_sym) << std::endl;
            std::cout << "error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}


#endif // _TEST_LIE_SPECIAL_ORTHOGONAL_H_

