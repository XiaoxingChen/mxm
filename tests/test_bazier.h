#if !defined(_TEST_BAZIER_H)
#define _TEST_BAZIER_H

#include "test_config.h"
#if TEST_AVAILABLE_BAZIER

#include "mxm/bazier.h"

using namespace mxm;


inline void testBazier()
{
    {
        if(factorial(3) != 6)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    {
        Vec x({30,30,23,12.5,0});
        Vec t({0.,0.25,0.5,0.75,1.});
        size_t n = 4;

        Vec expected_fx({0.,11.6484375,21.125,27.5859375,30.});
        if((expected_fx - bazier::f(n, x, t)).norm() > eps())
        {
            std::cout << mxm::to_string(bazier::f(n, x, t).T()) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

#else
inline void testBazier(){}
#endif

#endif // _TEST_BAZIER_H
