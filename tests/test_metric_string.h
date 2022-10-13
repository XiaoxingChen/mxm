#if !defined(_TEST_METRIC_STRING_H_)
#define _TEST_METRIC_STRING_H_

#include "mxm/metric_string.h"
namespace mxm{

void testLongestCommonString01()
{
    std::vector<size_t> s1{1,2,5,2,3};
    std::vector<size_t> s2{2,2,3,4};

    auto path = longestCommonString(s1, s2);
    for(auto & idx_pair : path)
    {
        if(s1.at(idx_pair[0]) != s2.at(idx_pair[1]))
        {
            for(auto & pair: path)
            {
                std::cout << "(" << s1.at(pair[0]) << "," << s2.at(pair[1]) << ")";
                throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
            }
        }
    }
}

void testLongestCommonString02()
{
    std::string s1("ABCDGH");
    std::string s2("AEDFHR");

    auto path = longestCommonString(s1, s2);
    std::string common_subseq;
    for(auto & idx_pair : path)
    {
        common_subseq.push_back(s1.at(idx_pair[0]));
    }
    if(common_subseq != std::string("ADH"))
    {
        std::cout << "path: " << mxm::to_string(path) << std::endl;
        std::cout << common_subseq << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testMetricString()
{
    testLongestCommonString01();
    testLongestCommonString02();
}

} //scope of namespace mxm
#endif // _TEST_METRIC_STRING_H_
