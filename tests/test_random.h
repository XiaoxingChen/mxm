#if !defined(_TEST_RANDOM_H_)
#define _TEST_RANDOM_H_

#include "mxm/random.h"
using namespace mxm;

void testRandom()
{
    for(size_t dim = 2; dim < 5; dim++)
    {
        auto vec = random::unitSphere<FloatType>(dim);
        if(vec.shape(0) != dim)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

#endif // _TEST_RANDOM_H_
