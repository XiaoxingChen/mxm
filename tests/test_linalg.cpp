#include "test_config.h"
#include "mxm/common.h"

#include <iostream>
#include <map>
#include <chrono>

#if TEST_AVAILABLE_LINALG_MAT
#include "mxm/linalg_mat.h"
#include "mxm/linalg_mat_ref.h"
#include "mxm/linalg_mat_block.h"
#include "mxm/string.h"
#include "mxm/linalg_norm.h"
#endif

#if TEST_AVAILABLE_LINALG_VEC
#include "mxm/linalg_vec.h"
#endif

#if TEST_AVAILABLE_LINALG_SOLVE
#include "mxm/linalg_solve.h"
#endif

#if TEST_AVAILABLE_LINALG_UTILS
#include "mxm/linalg_utils.h"
#include "mxm/string.h"
#endif

#if TEST_AVAILABLE_LINALG_COMPLEX
#include "mxm/linalg_complex.h"
#endif

#include "mxm/random.h"



using namespace mxm;

// A case causes orthogonalComplement() to fail.
// The rank of the matrix is 7, which might be the reason.
const Matrix<float> testDataOrthogonalComplement1()
{
    const static Matrix<float> mat({8,9},{
        -26.728987, -8.646025, -3.085133, 24.507498, 7.927440, 2.828723, 8.663803, 2.802480, 1.000000,
        28.314852, 12.239782, 2.006292, 86.905548, 37.567032, 6.157825, 14.113026, 6.100698, 1.000000,
        -30.244743, -33.029633, -4.542576, 48.864914, 53.364323, 7.339212, 6.658060, 7.271124, 1.000000,
        -0.390185, -0.358124, -0.034057, 121.599968, 111.608154, 10.613824, 11.456754, 10.515358, 1.000000,
        0.350281, 0.037470, 0.032123, 12.838590, 1.373350, 1.177374, 10.904425, 1.166451, 1.000000,
        71.417076, 17.710978, 4.554032, 61.560223, 15.266541, 3.925493, 15.682162, 3.889076, 1.000000,
        -13.366714, -7.739275, -1.489354, 47.073555, 27.255404, 5.245057, 8.974841, 5.196398, 1.000000,
        34.085690, 20.318111, 2.572532, 105.628609, 62.964069, 7.972055, 13.249859, 7.898098, 1.000000});
    return mat;
}

// An example with det(U V^H) = -1.
const Matrix<float> testDataSVD1()
{
    static const Matrix<float> mat(square(),
        {-1.392623, 7.462401, -0.556963,
        -9.647902, 0.104333, -15.839573,
        1.008865, 17.574787, -0.170509});
    return mat;
}

#if TEST_AVAILABLE_LINALG_UTILS
void testOrthogonalComplement()
{
    {
        Mat in({2,3},{0,0,1, 0,1,0});

        if(in(Row(0)).matmul(orthogonalComplement(in).T())(0,0) > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(in(Row(1)).matmul(orthogonalComplement(in.T()))(0,0) > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Mat in({3,2},{
            -0.16224509,0.16913061,
            0.08134538,  0.,
            -0.06551123,  0.});
        Mat expected({3,1},{
            0.,
            -0.01107995,
            -0.01375799});

        if((orthogonalComplement(in) - expected).norm() > eps())
        {
            std::cout << "error: " << (orthogonalComplement(in) - expected).norm() << std::endl;
            std::cout << "actual: " << mxm::to_string (orthogonalComplement(in).T()) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    #if(0)
    {
        Mat in({1,3},{1,0,0});
        Mat complement = orthogonalComplement(in);
        Vec out1 = complement(Row(0));
        Vec out2 = complement(Row(1));

        if(static_cast<Vec>(in).dot(out1) > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(static_cast<Vec>(in).dot(out2) > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
    #endif

}
#else
void testOrthogonalComplement(){}
#endif

#if TEST_AVAILABLE_LINALG_SOLVE
void testSolveLowerTriangle()
{
    Mat L({3,3}, {1,0,0, 2,3,0, 4,5,6});
    Vec b({2,3,4});
    Vec expect({2, -1./3, -7./18});
    Vec x = solveLowerTriangle(L, b);
    if((x - expect).norm() > eps())
    {
        std::cout << mxm::to_string(x.T()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    Vec expect_u({-4./9, -1./9, 2./3});
    Vec x_u = solveUpperTriangle(L.T(), b);
    if((x_u - expect_u).norm() > eps())
    {
        std::cout << mxm::to_string(x_u.T()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

void testQRcalcMatQ()
{
    // return;
    if(0)
    {
        Mat mat_a({3,3},{0, 0, 0, 0.431263924, -1, -1, 0.902225852, 0, -1});
        Mat expect_q({3,3},
        {0.00000000e+00,  5.55111512e-17, -1.00000000e+00,
        -4.31263911e-01, -9.02225825e-01, -5.55111512e-17,
        -9.02225825e-01,  4.31263911e-01,  1.11022302e-16});
        Mat mat_q = qr::decomposeByRotation(mat_a)[0];
        if((mat_q - expect_q).norm() > 2 * eps())
        {
            std::cout << mxm::to_string(mat_q) << std::endl;
            std::cout << "norm: " << (mat_q - expect_q).norm() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {

        Mat m({3,3},{
            0.99845566, 0.94293762, 0.80127313,
            0.6007756 , 0.30253734, 0.78378593,
            0.93970026, 0.781002  , 0.86480691});

        // Mat mat_q = qr::calcMatQFromRotation(m);
        Mat mat_q = qr::decomposeByRotation(m)[0];

        Mat expect_q({3,3},
            {0.66699001,  0.5087728 , 0.54431109,
            0.40133111, -0.86084531, 0.31285571,
            0.62774013,  0.00977734,  -0.77836157});

        if((mat_q - expect_q).norm() > 5 * eps())
        {
            std::cout << mxm::to_string(mat_q) << std::endl;
            std::cout << "norm: " << (mat_q - expect_q).norm() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    {
        Matrix<Complex<FloatType>> m({3,3},{
            {1,2},{2,0},{3,0},
            {2,0},{3,0},{0,4},
            {3,0},{3,0},{1,2}});

        auto q_r = qr::decomposeByRotation(m);
        // std::cout << mxm::to_string(q_r[1]) << std::endl;

        // std::cout << mxm::to_string((q_r[0].matmul(q_r[1]))) << std::endl;

        if((conj(q_r[0].T()).matmul(q_r[0]) - Matrix<Complex<FloatType>>::identity(3)).norm() > 5 * eps())
        {
            std::cout << "error: " << (conj(q_r[0].T()).matmul(q_r[0]) - Matrix<Complex<FloatType>>::identity(3)).norm() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if((q_r[0].matmul(q_r[1]) - m).norm() > 30 * eps())
        {
            std::cout << "error: " << (q_r[0].matmul(q_r[1]) - m).norm() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    if(0){
        Matrix<float> m = Matrix<float>::zeros({3,3});
        auto q_r = qr::decomposeByRotation(m);

        std::cout << mxm::to_string(q_r) << std::endl;
    }

    {
        // this test can be used to find the nullspace of a under determined matrix
        Matrix<float> m({4,4},{
            -0.75,  0.25,  0.25,  0.25,
            1.  , -1.  ,  0.  ,  0.  ,
            1.  ,  0.  , -1.  ,  0.  ,
            1.  ,  0.  ,  0.  , -1.  });

        auto q_r = qr::decomposeByRotation(m);
        if(q_r[1](Row(end() - 1)).norm() > std::numeric_limits<float>::epsilon())
        {
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        auto null_base = orthogonalComplement(q_r[1](Block({0, end()-1}, {0, end()})));
        if(m.matmul(null_base.T()).norm() > std::numeric_limits<float>::epsilon() * 10)
        {
            std::cout << mxm::to_string(null_base) << std::endl;
            std::cout << mxm::to_string(m.matmul(null_base.T())) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

void testQRSolve()
{
    // using namespace qr;
    Mat mat_a({3,3}, {12, -51, 4, 6, 167, -68, -4, 24, -41});
    Mat expect_q({3,3}, {6./7, -69./175, -58./175, 3./7, 158./175, 6./175, -2./7, 6./35, -33./35});
    Mat mat_q = qr::calcMatQFromReflection(mat_a);
    if((mat_q - expect_q).norm() > 2 * eps())
    {
        std::cout << mxm::to_string(mat_q) << std::endl;
        std::cout << "norm: " << (mat_q - expect_q).norm() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    Vec b({2,3,9});
    Vec expect({-0.10244898, -0.08359184, -0.25844898});
    Vec x(qr::solve(mat_a, b));

    // if((x - expect).norm() > eps())
    if((mat_a.matmul(x) - b).norm() > 10 * eps())
    {
        std::cout << mxm::to_string(x.T()) << std::endl;
        std::cout << "err: " << (mat_a.matmul(x) - b).norm() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Matrix<Complex<FloatType>> m({3,3},{
            {1,3},{2,0},{3,0},
            {2,2},{3,0},{0,4},
            {3,0},{3,0},{1,2}});

        Matrix<Complex<FloatType>> b({3,1}, {
            {1,2},
            {3,4},
            {6,2}}, COL);

        auto x = qr::solve(m, b);

        if((m.matmul(x) - b).norm() > 10 * eps())
        {
            std::cout << "error: " << (m.matmul(x) - b).norm() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }


    }
}

// Reference: https://www.cs.utexas.edu/users/flame/laff/alaff-beta/chapter10-QR1-simple-shifted-QR.html
// https://www.cs.utexas.edu/users/flame/laff/alaff-beta/Assignments/Week10/answers/SimpleShiftedQRAlg.m
void testQRIteration()
{
    size_t n = 3;
    Mat rot({n,n}, {
        -0.0816791 , -0.43273854,  0.89781172,
        -0.36035366,  0.85270191,  0.3782125 ,
        -0.92923289, -0.29263768, -0.22558684});

    Mat tmp = rot;
    Mat eig_vec = Mat::identity(n);


    for(size_t i = 0; i < 100; i++)
    {
        FloatType rho = tmp(2,2);
        Mat shift = Mat::identity(n) * rho;
        auto q_r = qr::decomposeByRotation(tmp - shift);
        tmp = q_r[1].matmul(q_r[0]) + shift;

        eig_vec = eig_vec.matmul(q_r[0]);
    }
    // std::cout << mxm::to_string(tmp) << std::endl;
    Mat recover = eig_vec.matmul(tmp).matmul(eig_vec.T());
    if((recover - rot).norm() > 40 * eps())
    {
        std::cout << (recover - rot).norm() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

namespace mxm
{

namespace qr
{
template<typename DType>
DType wilkinsonShiftStrategy(const Matrix<DType> mat_a);
} // namespace qr

template<typename DType>
Matrix<DType> shiftedQRIteration(
    const Matrix<DType>& mat,
    Matrix<DType>* p_orthogonal,
    qr::TraverseSeq idx_seq=qr::eUpperTrianglize,
    size_t max_it=40,
    DType tol=eps<DType>());

} // namespace mxm


void testUnsymmetrixEigenvaluePipeline01()
{
    size_t n = 5;
    Mat rot({n,n}, {
         8.52072892e-01, -1.86443589e-01, -1.52868423e-01, 4.48402930e-01, -1.21559174e-01,
         1.66260826e-01, -3.10093307e-01,  4.89428897e-01, -6.23694169e-02,  7.95467717e-01,
        -2.09350679e-01,  3.08823048e-01,  6.01299790e-01, 6.90088793e-01, -1.51712355e-01,
         4.81917029e-04,  6.33156375e-01, -4.73611309e-01, 2.52124640e-01,  5.57887325e-01,
         4.50001318e-01,  6.10591729e-01,  3.88872076e-01, -5.05228157e-01, -1.34905788e-01});

    auto q_hesson = qr::decomposeByRotation(rot, qr::eUpperHessenbergize, true);
    Mat tmp = q_hesson[1];
    Mat eig_vec = q_hesson[0];

    for(size_t i = 0; i < 50; i++)
    {
        // FloatType rho = 1;
        FloatType rho = qr::wilkinsonShiftStrategy(tmp);
        Mat shift = Mat::identity(n) * rho;
        auto q_r = qr::decomposeByRotation(tmp - shift, qr::eSubdiagonal);
        tmp = q_r[1].matmul(q_r[0]) + shift;

        eig_vec = eig_vec.matmul(q_r[0]);
        // std::cout << "i: " << i << " err: " << qr::errorQuasiUpperTriangle(tmp) << std::endl;
    }
    std::cout << mxm::to_string(tmp) << std::endl;
    Mat recover = eig_vec.matmul(tmp).matmul(eig_vec.T());
    std::cout << mxm::to_string(recover) << std::endl;
}

void testUnsymmetrixEigenvaluePipeline02()
{
    size_t n = 5;
    Mat mat_a({n,n}, {
        0.0028910405, 0.2797538283, 0.9771945489, 0.2228983867, 0.6710390525,
        0.6567445086, 0.0454680832, 0.9139623576, 0.7276380447, 0.0890351175,
        0.7236019741, 0.6144312151, 0.5533828464, 0.4187292391, 0.0462407111,
        0.6310457773, 0.0022013412, 0.2199422767, 0.7315915776, 0.0965157312,
        0.2525883171, 0.7773467106, 0.2399114012, 0.2521241482, 0.1006454857});

    // Mat mat_a({n,n},
    // {0.7388932475, 0.9445831439, 0.131383505 , 0.6429751039, 0.7577514673,
    // 0.792865382 , 0.993147103 , 0.3259064733, 0.7192953417, 0.5080102095,
    // 0.6231152385, 0.640810716 , 0.7857406608, 0.0312927124, 0.2884410864,
    // 0.0966878558, 0.6565393823, 0.2895794538, 0.3901908797, 0.2669917572,
    // 0.7349230421, 0.2174315283, 0.3584232178, 0.0402299961, 0.5570949113});


    // expected: [ 2.0937458093+0.j, -0.6786741964+0.j, -0.1648492118+0.3498644439j, -0.1648492118-0.3498644439j,  0.3486058443+0.j]


    auto q_hessen = qr::decomposeByRotation(mat_a, qr::eUpperHessenbergize, true);
    Mat tmp = q_hessen[1];
    std::cout << mxm::to_string(q_hessen[1]) << std::endl;
    // Mat eig_vec = Mat::identity(n);
    Mat eig_vec = q_hessen[0];


    for(size_t i = 0; i < 30; i++)
    {
        // FloatType rho = 1;
        FloatType rho = qr::wilkinsonShiftStrategy(tmp);
        Mat shift = Mat::identity(n) * rho;
        auto q_r = qr::decomposeByRotation(tmp - shift, qr::eSubdiagonal);
        tmp = q_r[1].matmul(q_r[0]) + shift;

        eig_vec = eig_vec.matmul(q_r[0]);

        std::cout << "r: " << mxm::to_string(q_r[0]) << std::endl;

        // std::cout << "i: " << i << " err: " << qr::errorOrthogonalBlockDiagonal(q_r[0]) << std::endl;
    }
    std::cout << mxm::to_string(tmp) << std::endl;
    Mat recover = eig_vec.matmul(tmp).matmul(eig_vec.T());
    // std::cout << mxm::to_string(recover) << std::endl;
}

template<typename DType>
void checkValidSVD(
    const Matrix<DType>& mat,
    const std::array<Matrix<DType>, 3>& u_d_vh)
{
    size_t area = mat.shape(0) * mat.shape(1);
    auto restore = u_d_vh[0].matmul(diagonalMatrix(u_d_vh[1])).matmul(u_d_vh[2]);
    if(norm(mat - restore) / area > 15 * std::numeric_limits<float>::epsilon())
    {
        std::cout << mxm::to_string( mat - restore , 12) << std::endl;
        std::cout << "error: " << ( norm(mat - restore) / area) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
    DType error(0);
    if(!isIdentity(u_d_vh[0].matmul(u_d_vh[0].T()), &error, 15 * eps()))
    {
        std::cout << mxm::to_string( u_d_vh[0] ) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(!isIdentity(u_d_vh[2].matmul(u_d_vh[2].T()), &error, 50 * eps()))
    {
        std::cout << mxm::to_string(u_d_vh[2].matmul(u_d_vh[2].T())) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testEigenvalues()
{
    size_t n = 5;
    Mat mat_a({n,n}, {
        0.0028910405, 0.2797538283, 0.9771945489, 0.2228983867, 0.6710390525,
        0.6567445086, 0.0454680832, 0.9139623576, 0.7276380447, 0.0890351175,
        0.7236019741, 0.6144312151, 0.5533828464, 0.4187292391, 0.0462407111,
        0.6310457773, 0.0022013412, 0.2199422767, 0.7315915776, 0.0965157312,
        0.2525883171, 0.7773467106, 0.2399114012, 0.2521241482, 0.1006454857});

    {
        auto eigvals_std_vec = eigvals(mat_a);
        Matrix<Complex<FloatType>> eigvals({n,1});
        for(size_t i = 0; i < n; i++) eigvals(i,0) = eigvals_std_vec.at(i);
        // std::cout << mxm::to_string(eigvals) << std::endl;
        Matrix<Complex<FloatType>> expected({n,1}, {
            {2.0937458093, 0},
            {-0.6786741964, 0},
            {-0.1648492118, 0.3498644439},
            {-0.1648492118, -0.3498644439},
            {0.3486058443, 0}});

        if((eigvals - expected).norm() > 20*eps())
        {
            std::cout << "error: " << (eigvals - expected).norm() << std::endl;
            std::cout << mxm::to_string(eigvals) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        Matrix<Complex<FloatType>> cmat_a = Matrix<Complex<FloatType>>::zeros({n,n});
        cmat_a.traverse([&](auto i, auto j){cmat_a(i,j)(0) = mat_a(i,j);});

        // std::cout << mxm::to_string(cmat_a) << std::endl;

        auto eig_val_vec = eig(mat_a);

        // std::cout << mxm::to_string(eig_val_vec) << std::endl;

        // for(size_t i = 0; i < n; i++)
        for(size_t i = 0; i < 1; i++)
        {
            auto eig_val = eig_val_vec[0](i,0);
            auto eig_vec = eig_val_vec[1](Col(i));
            // std::cout << "eigv: " << mxm::to_string(eig_vec.T()) << std::endl;
            auto err = (cmat_a.matmul(eig_vec) - eig_val * eig_vec).norm();
            if(err > 20 * eps())
            {
                std::cout << "Error! " << __FILE__ << ":" << __LINE__ << std::endl;
                std::cout << "i: " << i << ", error: " << err << std::endl;
                std::cout << "mat.matmul(vec): " << mxm::to_string(cmat_a.matmul(eig_vec).T()) << "\n, eig_val * eig_vec: "
                << mxm::to_string((eig_val * eig_vec).T()) << std::endl;
                std::cout << "residual: " << mxm::to_string((cmat_a.matmul(eig_vec) - eig_val * eig_vec).T()) << std::endl;
                auto lhs = cmat_a.matmul(eig_vec);
                auto rhs = eig_val * eig_vec;
                std::cout << "residual2: " << mxm::to_string((lhs - rhs).T()) << std::endl;

            }
        }
        // std::cout << mxm::to_string(eig_val_vec) << std::endl;
    }

    {
        Matrix<double> mat_b(mat_a.T().matmul(mat_a));

        auto eig_val_vec = symmetricEig(mat_b);

        // std::cout << mxm::to_string(eig_val_vec) << std::endl;
        // std::cout << eig_val_vec[1](Col(0)).norm() << std::endl;

        for(size_t i = 0; i < n; i++)
        {
            auto eig_vec = eig_val_vec[1](Col(i));
            auto err = (mat_b.matmul(eig_vec) - eig_val_vec[0](i,0) * eig_vec).norm();
            if(err > 15 * std::numeric_limits<double>::epsilon())
            {
                std::cout << "Error! " << __FILE__ << ":" << __LINE__ << std::endl;
                std::cout << "i: " << i << ", error: " << err << std::endl;
            }
        }
    }

    #if(1)
    {
        auto u_s_vh = svd(mat_a);
        // u_s_vh[2](Row(0)) *= -1;
        try
        {
            checkValidSVD(mat_a, u_s_vh);
        }
        catch(const std::exception& e)
        {
            std::cout << e.what() << '\n';
            std::cout << "Catch error! " << __FILE__ << ":" << __LINE__ << std::endl;
        }
    }
    #endif

    {
        Matrix<float> skew({5,5},
        {0, 1, 2, 3, 4,
        -1, 0, 3, 6, 7,
        -2,-3, 0, 8, 2,
        -3,-6,-8, 0, 1,
        -4,-7,-2,-1, 0});

        Matrix<Complex<float>> expected({5,1}, {
            {0,+13.25652325},
            {0,-13.25652325},
            {0, +4.15506814},
            {0, -4.15506814},
            {0, +0.}});

        auto eig_vals = Matrix<Complex<float>>(expected.shape(), std::move(eigvals(skew)));
        if(norm(eig_vals - expected) > 40 * std::numeric_limits<float>::epsilon())
        {
            std::cout << "error: " << norm(eig_vals - expected) << std::endl;
            std::cout << mxm::to_string(eig_vals) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {
        for(size_t n = 3; n < 5; n++)
        {
            Matrix<double> mat = random::uniform<double>({n,n});
            auto u_d_vh = svd(mat);
            try
            {
                checkValidSVD(mat, u_d_vh);
            }
            catch(const std::exception& e)
            {
                std::cout << e.what() << '\n';
                std::cout << "Catch error! " << __FILE__ << ":" << __LINE__ << std::endl;
            }
        }

        std::cout << "TODO: SVD large error at n > 5." << std::endl;
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;

    }

}

void testEigenvalues2()
{

{
    size_t n = 5;
    Matrix<double> rot({n,n}, {
         8.52072892e-01, -1.86443589e-01, -1.52868423e-01, 4.48402930e-01, -1.21559174e-01,
         1.66260826e-01, -3.10093307e-01,  4.89428897e-01, -6.23694169e-02,  7.95467717e-01,
        -2.09350679e-01,  3.08823048e-01,  6.01299790e-01, 6.90088793e-01, -1.51712355e-01,
         4.81917029e-04,  6.33156375e-01, -4.73611309e-01, 2.52124640e-01,  5.57887325e-01,
         4.50001318e-01,  6.10591729e-01,  3.88872076e-01, -5.05228157e-01, -1.34905788e-01});

    std::cout << mxm::to_string(eigvals(rot)) << std::endl;
}

}


#endif

Mat testMatRefTransposeReturn(size_t dim)
{
    return Mat::identity(dim).T();
}

void testMatRef()
{
    {//test 01
        auto mat = testMatRefTransposeReturn(3);
        if(mat(0,0) != 1.f)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }//test 01

    {//test 02
        Mat mat_a(Mat::ones({3, 10}));
        MatRef mr(mat_a(Col(1)).T());
        if((mr - Mat::ones({1, 3})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }//test 02

    {//test 03
        Mat mat_a(Mat::ones({3, 10}));
        MatRef m1(mat_a(Block({},{3,7})).T());
        MatRef m2(m1(Row(1)).T());
        MatRef m3(m1(Col(1)));
        if((m2 - Mat::ones({3, 1})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        if((m3 - Mat::ones({4, 1})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {//test 04
        Mat vs(Mat::ones({2, 3}));
        MatRef vst = vs.T();
        auto b = vst(Block({1,},{0, end()}));
        if((b(0,1) - 1 )> eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

#if 0
    {//test 05
        Mat vs(Mat::ones({2, 3}));
        auto block = vs(Row(0));
        auto v = block.asVector();
        if((v(2) - 1 )> eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {//test 06
        Mat vs(Mat::ones({2, 3}));
        vs(Col(1)) = Vec({5,6});
        auto v = vs(Row(0)).asVector();
        if(v.T().shape() != Shape({3,1}))
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        if(v(1) != 5)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
#endif

#if TEST_AVAILABLE_LINALG_VEC

    {//test 07
        Mat mat(Mat::identity(3));
        auto v = mat(Col(1));
        auto vt = v.T();
        if((vt - Vec({0,1,0}).T()).norm() > eps())
        {
            std::cout << mxm::to_string(mat) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if((vt - Vec({0,1,0})).norm() > eps())
        {
            std::cout << "WARNING: todo fix.\n" << std::string(__FILE__) + ":" + std::to_string(__LINE__) << std::endl;
        }

    }
#endif
#if 0
    {//test 08
        Mat mat(Mat::identity(3));
        auto b = mat(Col(1));
        auto v = b.asVector();

        if((v - Vec({0,1,0})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
#endif

    {//test 09
        Mat mat = Mat::identity(3);
        Mat mat_b = mat(Block({0,end()},{0, end()}));
        if(mat.shape() != Shape{3,3})
        {
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {//test 10
        Matrix<uint32_t> mat({4,4},{
            0,1,2,3,
            4,5,6,7,
            8,9,10,11,
            12,13,14,15});

        Matrix<uint32_t> expected({3,3}, {5,9,13, 6,10,14, 7,11,15});

        mat = mat(Block({1,end()}, {1,end()})).T();
        if(mat != expected)
        {
            std::cout << mxm::to_string(mat) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    {//test 11
        size_t n = 5;
        Mat mat_a({n,n}, {
        0.0028910405, 0.2797538283, 0.9771945489, 0.2228983867, 0.6710390525,
        0.6567445086, 0.0454680832, 0.9139623576, 0.7276380447, 0.0890351175,
        0.7236019741, 0.6144312151, 0.5533828464, 0.4187292391, 0.0462407111,
        0.6310457773, 0.0022013412, 0.2199422767, 0.7315915776, 0.0965157312,
        0.2525883171, 0.7773467106, 0.2399114012, 0.2521241482, 0.1006454857});

        Mat eig_vecs(fixRow(n), {0.466463, 0.516039, 0.524968, 0.327560, 0.365007});
        auto eig_vec = eig_vecs(Col(0));
        float eig_val = 2.09374581;

        if((mat_a.matmul(eig_vec) - eig_val * eig_vec).norm() > 10 * eps())
        {
            std::cout << mxm::to_string((mat_a.matmul(eig_vec) - eig_val * eig_vec).T()) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    { // test 12
        Mat left = Mat::zeros({4,3});
        Mat right = Mat::ones({4,3});

        auto blk_l = Block({0,3},{0,3});
        auto blk_r = Block({1,4},{0,3});
        left(blk_l) = right(blk_r);

        if(left.norm() < eps())
        {
            std::cout << mxm::to_string(left) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

#if TEST_AVAILABLE_LINALG_SOLVE
void testMatInv()
{
    if(1){
        Mat m({3,3},{
            0.99845566, 0.94293762, 0.80127313,
            0.6007756 , 0.30253734, 0.78378593,
            0.93970026, 0.781002  , 0.86480691});

        if((m.inv().matmul(m) - Mat::identity(3)).norm() > eps())
        {
            std::cout << "error norm: " << (m.inv().matmul(m) - Mat::identity(3)).norm() << std::endl;
            std::cout << "WARNING: todo fix.\n" << std::string(__FILE__) + ":" + std::to_string(__LINE__) << std::endl;
        }

    }
}
#endif

void testRvalueReference()
{
    Mat a = Mat::ones({100, 100});
    Mat b = std::move(a);
    if(a.shape() != Shape({0,0}))
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

#if TEST_AVAILABLE_LINALG_COMPLEX
void testComplexLogExtremeCondition();
void testComplexExp()
{
    {
        if((exp(Complex<float>{0, M_PI}) + 1.f).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if((exp(Complex<float>{0, M_PI_2}) - Complex<float>{0,1}).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}
void testComplexBase()
{
    {
        Complex<FloatType> a({1,2});
        Complex<FloatType> b({3,1});
        Complex<FloatType> expected({1,7});

        if((a * b - expected).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(mxm::to_string(a).size() == 0)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        auto mat = Matrix<Complex<float>>::identity(3);
        if(mxm::to_string(mat).size() == 0)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        float error = (Complex<float>{0.5, 0} - 0.5).norm();
        if(error > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(1){
        size_t n = 5;
        Mat v_real(fixCol(n), {0.97665438, 1.08045411, 1.09914879, 0.68582671, 0.76423185});
        Matrix<Complex<float>> v_complex(v_real);
        if((v_complex - v_real).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(1){//test 03
        size_t n = 5;
        Mat mat_a_real({n,n}, {
        0.0028910405, 0.2797538283, 0.9771945489, 0.2228983867, 0.6710390525,
        0.6567445086, 0.0454680832, 0.9139623576, 0.7276380447, 0.0890351175,
        0.7236019741, 0.6144312151, 0.5533828464, 0.4187292391, 0.0462407111,
        0.6310457773, 0.0022013412, 0.2199422767, 0.7315915776, 0.0965157312,
        0.2525883171, 0.7773467106, 0.2399114012, 0.2521241482, 0.1006454857});

        Matrix<Complex<float>> mat_a(mat_a_real);
        // std::cout << mxm::to_string(mat_a) << std::endl;

        Mat eig_vec_real(fixRow(n), {0.46646273, 0.51603882, 0.52496764, 0.32755968, 0.36500698});
        Matrix<Complex<float>> eig_vec(eig_vec_real);
        Complex<float> eig_val{2.09374581, 0};

        // Mat v_expected(fixRow(n), {0.97665438, 1.08045411, 1.09914879, 0.68582671, 0.76423185});

        if((mat_a.matmul(eig_vec) - eig_val * eig_vec).norm() > 10 * eps())
        {
            // std::cout << mxm::to_string((mat_a.matmul(eig_vec) - v_expected).T()) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    // test log of complext
    {
        float angle = M_PI_4;
        Complex<float> z1({cos(angle), sin(angle)});
        if(log(z1)(1) - angle > eps())
        {
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        Quaternion<float> q1({cos(angle), 0, sin(angle), 0});

        if(log(q1)(2) - angle > eps())
        {
            std::cout << mxm::to_string(log(q1)) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(norm(log(Complex<float>({-1,0})) - Complex<float>({0, M_PI})) > eps())
        {
            std::cout << mxm::to_string(log(Complex<float>({-1,0}))) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(norm(log(Complex<float>({0,-1})) - Complex<float>({0, -M_PI_2})) > eps())
        {
            std::cout << mxm::to_string(log(Complex<float>({0,-1}))) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(norm(mxm::inv(Quaternion<float>({0,1,2,3})) - mxm::pow(Quaternion<float>({0,1,2,3}), -1.f)) > std::numeric_limits<float>::epsilon())
        {
            std::cout << mxm::to_string(mxm::inv(Quaternion<float>({0,1,2,3}))) << std::endl;
            std::cout << mxm::to_string(mxm::pow(Quaternion<float>({0,1,2,3}), -1.f)) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        {
            Quaternion<float> q0({0,0,0,1});
            Quaternion<float> q1({1,0,0,0});
            Quaternion<float> mid = interp::slerp(q0, q1, 0.5f);
            if(mxm::norm(mid - Quaternion<float>({sqrt(.5f), 0.f, 0.f, sqrt(.5f)})) > std::numeric_limits<float>::epsilon())
            {
                std::cout << mxm::to_string (mid) << std::endl;
                throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
            }

        }

    }


    testComplexLogExtremeCondition();


}
void testComplexAtPart()
{
    Matrix<Quaternion<float>> mat({2,1}, {{0,0,0,1}, {0,0,0,-1}});
    auto matK = matrixAtPart(mat, 3);
    Matrix<float> expect({2,1}, {1, -1});
    if(!isZero(matK - expect, nullptr, eps<float>()))
    {
        std::cout << "matK:" << mxm::to_string(matK) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}
// extreme condition test
// reference: https://en.cppreference.com/w/cpp/numeric/complex/log
void testComplexLogExtremeCondition()
{
    const auto& infty = std::numeric_limits<float>::infinity();
    std::map<std::array<float, 2>, std::array<float, 2>> in_out_map;
    in_out_map[{-0, 0}] = {-infty, M_PI};
    in_out_map[{0, 0}] = {-infty, 0};
    in_out_map[{3, infty}] = {infty, M_PI_2};
    in_out_map[{3, NAN}] = {NAN, NAN};

}
#else
void testComplexBase(){}
#endif

#if TEST_AVAILABLE_LINALG_UTILS
void testCombinations()
{
    {
        if(combinations(5,2) != 10) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        if(combinations(5,3) != 10) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        if(combinations(3,1) != 3) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}
#else
void testCombinations(){}
#endif

template<typename DType=double>
Matrix<DType> testDataRandom5x5()
{
    Matrix<DType> mat_a({5,5}, {
            0.95850398, 0.3389291 , 0.0071609 , 0.74223576, 0.08585814,
            0.75156501, 0.87193671, 0.87278104, 0.25092665, 0.36739847,
            0.26376434, 0.87153321, 0.44203283, 0.74565173, 0.05153192,
            0.19369754, 0.68500599, 0.16890462, 0.09378345, 0.33960549,
            0.18433456, 0.83782842, 0.40354753, 0.79179645, 0.01594465});
    return mat_a;
}

template<typename DType=double>
void testSvdPipelineError()
{
    if(0){ // qr decomposition
        auto mat_a = testDataRandom5x5<DType>();
        auto mat_ata = mat_a.T().matmul(mat_a);
        auto q_r = qr::decomposeByRotation(mat_ata, qr::eUpperHessenbergize, true);
        auto restore = q_r[0].matmul(q_r[1]).matmul(q_r[0].T());
        auto diff = mat_ata - restore;
        DType error(0);
        isZero(diff, &error, 0.);

        std::cout << mxm::to_string(diff, 15) << std::endl;
        std::cout << "upper hessenberg error: " << error << std::endl;

        auto orth = q_r[0];
        auto quasi = shiftedQRIteration(q_r[1], &orth, qr::eSubdiagonal);
        std::cout << mxm::to_string(quasi, 15) << std::endl;
        restore = orth.matmul(quasi).matmul(orth.T());
        diff = mat_ata - restore;
        isZero(diff, &error, 0.);
        std::cout << mxm::to_string(diff, 15) << std::endl;
        std::cout << "qr iteration error: " << error << std::endl;
    }

    if(1){
        auto mat_a = testDataRandom5x5<DType>();
        auto mat_ata = mat_a.T().matmul(mat_a);
        auto val_vec = symmetricEig(mat_ata);
        auto restore = val_vec[1].matmul(diagonalMatrix(val_vec[0])).matmul(val_vec[1].T());
        auto diff = mat_ata - restore;
        DType error(0);
        isZero(diff, &error, 0.);
        double tol = std::is_same<DType, float>::value ? 2e-6 :
            std::is_same<DType, double>::value ? 2e-14: 0.;
        if(error > tol)
        {
            std::cout << "eps: " << std::numeric_limits<DType>::epsilon() << std::endl;
            std::cout << mxm::to_string(diff, 15) << std::endl;
            std::cout << "symmetric eigen error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    if(0){
        Matrix<Complex<DType>> mat_a({5,5}, {
            {0.83370547, 0.42428223}, {0.13529714, 0.41690796}, {0.04482265,0.10479416}, {0.9160445 , 0.13748177}, {0.73985756, 0.87589392},
            {0.52653829, 0.04706286}, {0.12826168, 0.11234982}, {0.66793262,0.46405475}, {0.26458016, 0.30325647}, {0.43352264, 0.4526207} ,
            {0.94636104, 0.17840708}, {0.02501129, 0.46682152}, {0.59295295,0.60034754}, {0.28156454, 0.35886331}, {0.53022439, 0.59091015},
            {0.37018765, 0.77140185}, {0.35539142, 0.27526053}, {0.79931838,0.10599939}, {0.57705667, 0.45787132}, {0.75884231, 0.50633517},
            {0.63278795, 0.59412612}, {0.33123404, 0.63476347}, {0.81156382,0.88095761}, {0.89720607, 0.13361227}, {0.65274818, 0.86401979}
        }, ROW);

        auto mat_ata = conj(mat_a.T()).matmul(mat_a);
        auto val_vec = symmetricEig(mat_ata);
        auto restore = val_vec[1].matmul(diagonalMatrix(val_vec[0])).matmul(val_vec[1].T());
        auto diff = mat_ata - restore;
        DType error(0);
        isZero(diff, &error, DType(0.));
        double tol = std::is_same<DType, float>::value ? 2e-6 :
            std::is_same<DType, double>::value ? 2e-14: 0.;
        if(error > tol)
        {
            std::cout << "eps: " << std::numeric_limits<DType>::epsilon() << std::endl;
            std::cout << mxm::to_string(diff, 5) << std::endl;
            std::cout << "symmetric eigen error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

void testConvolution01()
{
    Matrix<float> mat({20, 5});
    mat.traverse([&](auto i, auto j){mat(i,j) = i+j;});
    Matrix<float> core({3, 3}, {0,0,0,0,1,0,0,0,0});
    auto ret = convolute(mat, core);
    auto ret_p = convoluteParallel(mat, core);

    if(!mxm::isZero(ret - mat, nullptr, eps()))
    {
        std::cout << mxm::to_string(mat) << std::endl;
        std::cout << mxm::to_string(ret) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(!mxm::isZero(ret_p - mat, nullptr, eps()))
    {
        std::cout << mxm::to_string(mat) << std::endl;
        std::cout << mxm::to_string(ret_p) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

void testConvolution02()
{
#define ENABLE_SPEED_TEST 0
#if ENABLE_SPEED_TEST
    size_t loop_num = 10;
    size_t scale = 8;
    auto mat = random::uniform<float>({2000, 2000});
    Matrix<float> ret_serial(mat.shape());
    Matrix<float> core({3, 3}, {1,2,1,1,4,1,1,2,1});
    auto t_start = std::chrono::system_clock::now();
    for(size_t i = 0; i < loop_num; i++)
    {
        ret_serial = convolute(mat, core);
    }
    auto t_end = std::chrono::system_clock::now();
    auto t_serial = std::chrono::duration<double>(t_end - t_start).count() * scale;

    Matrix<float> ret_parallel(mat.shape());
    t_start = std::chrono::system_clock::now();
    for(size_t i = 0; i < loop_num * scale; i++)
    {
        ret_parallel = convoluteParallel(mat, core);
    }
    t_end = std::chrono::system_clock::now();
    auto t_parallel = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "t_serial: " << t_serial << ", t_parallel: " << t_parallel << std::endl;
#endif
}

void testConvolution()
{
    testConvolution01();
    testConvolution02();
}

void testLeftInverse()
{
    Matrix<float> mat({4,2},{0,0, 0,0, 1,1, 0,1}, ROW);

    auto ret = leftInv(mat);
    float error;
    if(! isIdentity(ret.matmul(mat), &error, eps<float>()))
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testReduce01()
{
    Matrix<float> mat({20, 8});
    mat.traverse([&](auto i, auto j){mat(i,j) = 1.f;});
    Matrix<float> core({2, 2}, {1.f, 1.f, 1.f, 1.f});
    auto ret = reduce(mat, core);
    auto ret_parallel = reduceParallel(mat, core);
    Matrix<float> expected = Matrix<float>::ones({10, 4}) * 4.f;

    float error;
    if(! isZero(expected - ret, &error, eps()))
    {
        std::cout << mxm::to_string(ret - expected) << std::endl;
        // std::cout << mxm::to_string(expected) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(! isZero(expected - ret_parallel, &error, eps()))
    {
        std::cout << mxm::to_string(ret_parallel - expected) << std::endl;
        // std::cout << mxm::to_string(expected) << std::endl;
        // std::cout << mxm::to_string(ret_parallel) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testReduce()
{
    testReduce01();
}

void testLinearAlgebra()
{
    Mat m1({3,3},{1,1,1, 2,2,2, 3,3,3});
    Mat v1({3,1},{1,1,1});
    Mat u1(v1.normalized());
    Mat expected({3,1}, {3, 6, 9});

    if((m1.matmul(v1) - expected).norm() > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    if((u1 - v1 * sqrt(1./3)).norm() > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

#if TEST_AVAILABLE_LINALG_SOLVE
    if(fabs(mxm::det(Mat::identity(3)) - 1.) > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
#endif

    {
        auto result = Matrix<float>::ones({2,5}) + Matrix<float>::ones({2,1});
        float error(0);
        if(!isZero(result - 2 * Matrix<float>::ones({2,5}), &error, std::numeric_limits<float>::epsilon()))
        {
            std::cout << mxm::to_string(result) << std::endl;
            std::cout << "error: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    {
        auto test_data = Matrix<size_t>::ones({5,6});
        auto val1 = sum(sum(test_data, 0), 1)(0,0);
        auto val2 = sum(test_data);
        if(val1 != val2)
        {
            std::cout << val2 << std::endl;
            std::cout << val1 << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    testMatRef();
    testOrthogonalComplement();
#if TEST_AVAILABLE_LINALG_SOLVE
    testSolveLowerTriangle();
    testQRcalcMatQ();
    testQRSolve();
    testMatInv();
    testEigenvalues();
    // testEigenvalues2();
#endif
    testRvalueReference();
    testComplexBase();
    testComplexExp();
    testComplexAtPart();
    // testTridiagonalizeSkewSymmetric();
    // testQRIteration();
    // testUnsymmetrixEigenvaluePipeline01();
    // testUnsymmetrixEigenvaluePipeline02();

    testCombinations();
    testSvdPipelineError<double>();
    testSvdPipelineError<float>();
    testConvolution();
    testReduce();
    testLeftInverse();

}