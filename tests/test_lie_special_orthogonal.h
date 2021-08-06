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

template<typename DType>
Matrix<DType> testDataPose()
{
    auto rot = so::exp<3>(so::wedge(Vector<double>{2,1,1}));
    auto pose = Matrix<DType>::identity(4);
    pose(se::rotBlk<3>()) = rot;
    pose(se::traBlk<3>()) = Vector<DType>{1,2,3};
    return pose;
}

template<typename DType>
Matrix<DType> testDataPose1()
{
    auto rot = so::exp<3>(so::wedge(Vector<double>{-1,2,1}));
    auto pose = Matrix<DType>::identity(4);
    pose(se::rotBlk<3>()) = rot;
    pose(se::traBlk<3>()) = Vector<DType>{1,1,0.5};
    return pose;
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

    {
        using DType = double;
        auto rot = so::exp<3>(so::wedge(Vector<double>{2,1,1}));
        auto pose = Matrix<DType>::identity(4);
        pose(se::rotBlk<3>()) = rot;
        pose(se::traBlk<3>()) = Vector<DType>{1,2,3};

        auto mid_pose = SE::interp(Matrix<DType>::identity(4), pose, 0.5);
        auto pose_sym = mid_pose.matmul(mid_pose);

        DType error(0);
        if(!isZero(pose - pose_sym, &error, 10*eps<DType>()))
        {
            std::cout << mxm::to_string(pose) << std::endl;
            std::cout << mxm::to_string(pose_sym) << std::endl;
            std::cout << "error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    { //test mat_q
        using DType = double;
        auto pose = testDataPose<DType>();
        const size_t LARGE_NUM = 20;
        auto mat_q_by_sum = Matrix<DType>::zeros({3,3});
        auto pose_alg = SE::log(pose);
        auto rho = so::wedge(pose_alg(se::traBlk()));
        auto phi = pose_alg(se::rotBlk());
        auto phi_n = Matrix<DType>::identity(3);
        for(size_t n = 0; n < LARGE_NUM; n++)
        {
            auto tmp = phi_n.matmul(rho) * (DType(1) / factorial(n+2));;
            for(size_t m = 0; m < LARGE_NUM; m++)
            {
                mat_q_by_sum += tmp;
                tmp = tmp.matmul(phi) * (DType(1) / (n+m+3));
            }
            phi_n = phi_n.matmul(phi);
        }

        auto mat_q = se::matQ(pose_alg);
        DType error(0);
        if(!isZero(mat_q_by_sum - mat_q, &error, 1e-10))
        {
            std::cout << "mat_q_by_sum:\n" << mxm::to_string(mat_q_by_sum) << std::endl;
            std::cout << "mat_q:\n" << mxm::to_string(mat_q) << std::endl;
            std::cout << "error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    { //test derive distance
        using DType = double;
        const size_t N = 3;
        auto pose1 = testDataPose<DType>();
        auto pose2 = testDataPose1<DType>();

        auto deriv = SE::derivDistance(pose1, SE::inv(pose2));

        DType delta = 1e-6;
        auto perturb = Matrix<DType>::identity(2*N) * delta;
        auto numerical_deriv = Matrix<DType>::identity(2*N);

        auto actual_cost = SE::log(pose1.matmul(SE::inv(pose2)));
        for(size_t i = 0; i < perturb.shape(1); i++)
        {
            auto perturbed_pose1 = se::exp(se::wedge(perturb(Col(i)))).matmul(pose1);
            auto perturbed_cost = SE::log(perturbed_pose1.matmul(SE::inv(pose2)));
            numerical_deriv(Col(i)) = se::vee(perturbed_cost - actual_cost) / delta;
        }

        DType error(0);
        if(!isZero(numerical_deriv - deriv, &error, 1e-5))
        {
            std::cout << "numerical_deriv:\n" << mxm::to_string(numerical_deriv) << std::endl;
            std::cout << "deriv:\n" << mxm::to_string(deriv) << std::endl;
            std::cout << "error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}


#endif // _TEST_LIE_SPECIAL_ORTHOGONAL_H_

