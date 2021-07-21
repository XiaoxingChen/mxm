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

inline void testLieSpecialUnitary()
{
    {
        auto mat = testDataSU2();
        float error(0);
        bool ret = SU::isValid(mat, std::numeric_limits<float>::epsilon(), &error);
        if(!ret)
        {
            std::cout << "is SU2 err: " << error << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
        
    }
}


#endif // _TEST_LIE_SPECIAL_UNITARY_H_
