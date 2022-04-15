#if !defined(_TEST_METRIC_STRING_H_)
#define _TEST_METRIC_STRING_H_

#include "mxm/metric_string.h"
using namespace mxm;

void testLongestCommonString01()
{
    std::vector<size_t> s1{1,2,5,2,3};
    std::vector<size_t> s2{2,2,3,4};

    size_t l = longestCommonString(s1, s2);
    if(l != 3)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testMetricString()
{
    testLongestCommonString01();
}


#endif // _TEST_METRIC_STRING_H_
