#if !defined(_TEST_LIE_SPECIAL_UNITARY_H_)
#define _TEST_LIE_SPECIAL_UNITARY_H_

#include "mxm/lie_special_unitary.h"

using namespace mxm;

inline Matrix<Complex<float>> testDataSU2()
{
    Vector<float> v{0.23570226, 0.47140452, 0.70710678, 0.47140452};
    Complex<float> a{v(0), v(1)};
    Complex<float> b{v(2), v(3)};
    return Matrix<Complex<float>>({2,2}, {a, -b.conj(), b, a.conj()}, ROW);
}

inline Matrix<Complex<float>> testDataSU2_01()
{
    Vector<float> v{0.5, 0.5, 0.5, 0.5};
    Complex<float> a{v(0), v(1)};
    Complex<float> b{v(2), v(3)};
    return Matrix<Complex<float>>({2,2}, {a, -b.conj(), b, a.conj()}, ROW);
}

inline void testLieSpecialUnitary()
{
    {
        auto mat = testDataSU2();
        float error(0);
        if(!SU::isValid(mat, &error, eps<float>()))
        {
            std::cout << "is SU2 err: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    {
        auto mat = testDataSU2();
        auto lie_alg = SU::log<2>(mat);
        auto mat_eq = su::exp<2>(lie_alg);
        float error;
        if(!isZero(mat - mat_eq, &error))
        {
            std::cout << "error: " << error << std::endl;
            std::cout << mxm::to_string(mat) << std::endl;
            std::cout << mxm::to_string(lie_alg) << std::endl;
            std::cout << mxm::to_string(mat_eq) << std::endl;
        }
    }
}


#endif // _TEST_LIE_SPECIAL_UNITARY_H_
