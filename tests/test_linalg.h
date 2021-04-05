#if !defined(_TEST_LINALG_H)
#define _TEST_LINALG_H

#include <iostream>
#include "mxm/linalg.h"
// #include "pixel.h"
// #include "accessor.h"


using namespace mxm;

inline void testOrthogonalComplement()
{
    {
        Mat in({2,3},{0,0,1, 0,1,0});

        if(in(Row(0)).matmul(orthogonalComplement(in).T())(0,0) > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(in(Row(1)).matmul(orthogonalComplement(in.T()))(0,0) > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Mat in({3,2},{-0.16224509,0.16913061, 0.08134538,  0., -0.06551123,  0.});
        Mat expected({3,1},{0., -0.01107995, -0.01375799});

        if((orthogonalComplement(in) - expected).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(0){
        Mat in({1,3},{1,0,0});
        Mat complement = orthogonalComplement(in);
        Vec out1 = complement(Row(0));
        Vec out2 = complement(Row(1));

        if(static_cast<Vec>(in).dot(out1) > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(static_cast<Vec>(in).dot(out2) > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testSolveLowerTriangle()
{
    Mat L({3,3}, {1,0,0, 2,3,0, 4,5,6});
    Vec b({2,3,4});
    Vec expect({2, -1./3, -7./18});
    Vec x = solveLowerTriangle(L, b);
    if((x - expect).norm() > eps())
    {
        std::cout << x.T().str() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    Vec expect_u({-4./9, -1./9, 2./3});
    Vec x_u = solveUpperTriangle(L.T(), b);
    if((x_u - expect_u).norm() > eps())
    {
        std::cout << x_u.T().str() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testQRcalcMatQ()
{
    // return;
    if(0)
    {
        Mat mat_a({3,3},{0, 0, 0, 0.431263924, -1, -1, 0.902225852, 0, -1});
        Mat expect_q({3,3},
        {0.00000000e+00,  5.55111512e-17, -1.00000000e+00,
        -4.31263911e-01, -9.02225825e-01, -5.55111512e-17,
        -9.02225825e-01,  4.31263911e-01,  1.11022302e-16});
        Mat mat_q = qr::calcMatQ(mat_a);
        if((mat_q - expect_q).norm() > 2 * eps())
        {
            std::cout << mat_q.str() << std::endl;
            std::cout << "norm: " << (mat_q - expect_q).norm() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {

        Mat m({3,3},{
            0.99845566, 0.94293762, 0.80127313,
            0.6007756 , 0.30253734, 0.78378593,
            0.93970026, 0.781002  , 0.86480691});

        Mat mat_q = qr::calcMatQFromRotation(m);

        Mat expect_q({3,3},
            {0.66699001,  0.5087728 , 0.54431109,
            0.40133111, -0.86084531, 0.31285571,
            0.62774013,  0.00977734,  -0.77836157});

        if((mat_q - expect_q).norm() > 5 * eps())
        {
            std::cout << mat_q.str() << std::endl;
            std::cout << "norm: " << (mat_q - expect_q).norm() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }
}

inline void testQRSolve()
{
    // using namespace qr;
    Mat mat_a({3,3}, {12, -51, 4, 6, 167, -68, -4, 24, -41});
    Mat expect_q({3,3}, {6./7, -69./175, -58./175, 3./7, 158./175, 6./175, -2./7, 6./35, -33./35});
    Mat mat_q = qr::calcMatQFromReflection(mat_a);
    if((mat_q - expect_q).norm() > 2 * eps())
    {
        std::cout << mat_q.str() << std::endl;
        std::cout << "norm: " << (mat_q - expect_q).norm() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    Vec b({2,3,9});
    Vec expect({-0.10244898, -0.08359184, -0.25844898});
    Vec x(qr::solve(mat_a, b));

    if((x - expect).norm() > eps())
    {
        std::cout << x.T().str() << std::endl;
        // std::cout << "norm: " << (mat_q - expect_q).norm() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline Mat testMatRefTransposeReturn(size_t dim)
{
    return Mat::Identity(dim).T();
}

inline void testMatRef()
{
    {
        auto mat = testMatRefTransposeReturn(3);
        if(mat(0,0) != 1.f)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Mat mat_a(Mat::ones({3, 10}));
        MatRef mr(mat_a(Col(1)).T());
        if((mr - Mat::ones({1, 3})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Mat mat_a(Mat::ones({3, 10}));
        MatRef m1(mat_a(Block({},{3,7})).T());
        MatRef m2(m1(Row(1)).T());
        MatRef m3(m1(Col(1)));
        if((m2 - Mat::ones({3, 1})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        if((m3 - Mat::ones({4, 1})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Mat vs(Mat::ones({2, 3}));
        MatRef vst = vs.T();
        auto b = vst(Block({1,},{0, end()}));
        if((b(0,1) - 1 )> eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(1){
        Mat vs(Mat::ones({2, 3}));
        auto block = vs(Row(0));
        auto v = block.asVector();
        if((v(2) - 1 )> eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Mat vs(Mat::ones({2, 3}));
        vs(Col(1)) = Vec({5,6});
        auto v = vs(Row(0)).asVector();
        if(v.T().shape() != Shape({3,1}))
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        if(v(1) != 5)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Mat mat(Mat::Identity(3));
        auto v = mat(Col(1));
        auto vt = v.T();
        if((vt - Vec({0,1,0})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Mat mat(Mat::Identity(3));
        auto b = mat(Col(1));
        auto v = b.asVector();

        if((v - Vec({0,1,0})).norm() > eps())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testMatInv()
{
    if(1){
        Mat m({3,3},{
            0.99845566, 0.94293762, 0.80127313,
            0.6007756 , 0.30253734, 0.78378593,
            0.93970026, 0.781002  , 0.86480691});

        Mat inv_expect({3,3},
            {-125.49635084,  -67.90837523,  177.82291256,
            77.68518512,   39.56953927, -107.84037325,
            66.20745954,   38.05430802,  -94.67626698});

        if((m.inv() - inv_expect).norm()> 0.01)
        {
            std::cout << m.inv().str() << std::endl;
            std::cout << "Norm: " << (m.inv() - inv_expect).norm() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if((m.inv() - inv_expect).norm() > eps())
        {
            std::cout << "error norm: " << (m.inv() - inv_expect).norm() << std::endl;
            std::cout << "WARNING: todo fix.\n" << std::string(__FILE__) + ":" + std::to_string(__LINE__) << std::endl;
        }

    }
}

inline void testRvalueReference()
{
    Mat a = Mat::ones({100, 100});
    Mat b = std::move(a);
    if(a.shape() != Shape({0,0}))
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}

// inline void testComplexBase()
// {
//     Complex a({1,2});
//     Complex b({3,1});
//     Complex expected({1,7});

//     std::cout << a.str() << std::endl;
//     if((a * b - expected).norm() > eps())
//         throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
// }

inline void testLinearAlgebra()
{
    Mat m1({3,3},{1,1,1, 2,2,2, 3,3,3});
    Mat v1({3,1},{1,1,1});
    Mat u1(v1.normalized());
    Mat expected({3,1}, {3, 6, 9});

    if((m1.matmul(v1) - expected).norm() > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    if((u1 - v1 * sqrt(1./3)).norm() > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    if(fabs(Mat::Identity(3).det() - 1.) > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));


    testMatRef();
    testMatInv();
    testOrthogonalComplement();
    testSolveLowerTriangle();
    testQRcalcMatQ();
    testQRSolve();
    testRvalueReference();
    // testComplexBase();


}


#endif // _TEST_LINALG_H
