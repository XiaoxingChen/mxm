#include "test_config.h"
#include "mxm/linalg_dual_number.h"

using namespace mxm;

template<typename DType>
DType testFunction01(DType x)
{
    return 3.f*x*x + 2.f*x;
}

void dualNumberTest01()
{
    float x = 2.f;
    float dfdx = testFunction01(DualNumber<float>{x, 1})(1);
    if(abs(dfdx - 14) > eps())
    {
        std::cout << "dfdx: " << dfdx << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

// void testDualNumberAtPart()
// {
//     Matrix<DualNumber<float>> mat({2,1}, {{0,1}, {0,-1}});
//     auto matK = matrixAtPart(mat, 3);
//     Matrix<float> expect({2,1}, {1, -1});
//     if(!isZero(matK - expect, nullptr, eps<float>()))
//     {
//         std::cout << "matK:" << mxm::to_string(matK) << std::endl;
//         throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
//     }
// }

void testDualNumberTypeConversion()
{
    Matrix<DualNumber<float>> mat({3,3});
    mat = Matrix<float>::identity(3);
    if(!isZero( matrixAtPart(mat, 0)  - Matrix<float>::identity(3)))
    {
        std::cout << "mat: " << mxm::to_string(mat) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

void testDualNumberTypeCast()
{
    float a = DualNumber<float>(3.);
    if(abs(a - 3.) > 0)
    {
        std::cout << "a: " << a << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testDualNumber()
{
    dualNumberTest01();
    testDualNumberTypeConversion();
    testDualNumberTypeCast();
}
