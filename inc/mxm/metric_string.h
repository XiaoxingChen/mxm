#if !defined(_METRIC_STRING_H_)
#define _METRIC_STRING_H_

#include "linalg.h"

namespace mxm
{

template<typename DType>
// std::vector<std::array<size_t, 2>>
size_t
longestCommonString(const std::vector<DType>& s1, const std::vector<DType>& s2)
{
    Matrix<size_t> dp = Matrix<size_t>::zeros({s1.size(), s2.size()});
    for(size_t i = 1; i < s1.size(); i++)
    {
        if(s1.at(i) == s2.at(0)) dp(i, 0) = 1;
        else dp(i, 0) = dp(i-1, 0);
    }
    for(size_t i = 1; i < s2.size(); i++)
    {
        if(s1.at(0) == s2.at(i)) dp(0, i) = 1;
        else dp(0, i) = dp(0, i-1);
    }
    for(size_t i = 1; i < s1.size(); i++)
    {
        for(size_t j = 1; j < s2.size(); j++)
        {
            size_t max_inherit = std::max(dp(i, j-1), dp(i-1, j));
            if(s1.at(i) == s2.at(j))
            {
                max_inherit = std::max(max_inherit, dp(i-1, j-1) + 1);
            }
            dp(i,j) = max_inherit;
        }
    }
    // back trace
#if 0
    std::vector<std::array<size_t, 2>> path;

    size_t i = dp.shape(0) - 1;
    size_t j = dp.shape(1) - 1;
    path.push_back({i,j});
    size_t curr_val = dp(i,j);
    while(i > 0 && j > 0)
    {
        if(curr_val - 1 == dp(i-1, j-1))
        {
            path.push_back({i-1, j-1});
            continue;
        }

    }

    return path;
#else
    return dp(s1.size() - 1, s2.size() - 1);
#endif
}

} // namespace mxm


#endif // _METRIC_STRING_H_
