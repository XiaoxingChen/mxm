#if !defined(_METRIC_STRING_H_)
#define _METRIC_STRING_H_

#include "linalg.h"

namespace mxm
{

const size_t S1_INC = 1;
const size_t S2_INC = 2;
const size_t BOTH_INC = 3;

template<typename StrType> size_t
longestCommonStringHelper_(
    const StrType& s1, const StrType& s2,
    size_t end1, size_t end2,
    Matrix<size_t>& dp_table, Matrix<size_t>& pred_table)
{
    if(dp_table(end1, end2) != std::numeric_limits<size_t>::max())
    {
        return dp_table(end1, end2);
    }

    size_t s1_inc_reward = longestCommonStringHelper_(s1, s2, end1-1, end2, dp_table, pred_table);
    size_t s2_inc_reward = longestCommonStringHelper_(s1, s2, end1, end2-1, dp_table, pred_table);
    size_t both_inc_reward = s1.at(end1 - 1) == s2.at(end2 - 1) ? longestCommonStringHelper_(s1, s2, end1-1, end2-1, dp_table, pred_table) + 1 : 0;
    size_t max_reward = std::max(both_inc_reward, std::max(s1_inc_reward, s2_inc_reward));
    dp_table(end1, end2) = max_reward;

    pred_table(end1, end2) =
        (max_reward == s1_inc_reward) ? S1_INC
        : (max_reward == s2_inc_reward) ? S2_INC
        : BOTH_INC;

    return dp_table(end1, end2);
}

template<typename StrType>
std::vector<std::array<size_t, 2>>
longestCommonString(const StrType& s1, const StrType& s2)
{
    Matrix<size_t> dp = Matrix<size_t>::ones({s1.size()+1, s2.size()+1}) * std::numeric_limits<size_t>::max();
    Matrix<size_t> best_pred = Matrix<size_t>::zeros({s1.size()+1, s2.size()+1});
    for(size_t i = 0; i < s1.size()+1; i++)  dp(i, 0) = 0;
    for(size_t j = 0; j < s2.size()+1; j++)  dp(0, j) = 0;

    longestCommonStringHelper_(s1, s2, s1.size(), s2.size(), dp, best_pred);
    // std::cout << mxm::to_string(best_pred) << std::endl;
    // std::cout << mxm::to_string(dp) << std::endl;
    // back trace

    std::vector<std::array<size_t, 2>> path;

    size_t i = dp.shape(0) - 1;
    size_t j = dp.shape(1) - 1;
    while(i > 0 && j > 0)
    {
        // path.push_back({i,j});
        if(best_pred(i,j) == S1_INC)
        {
            i--;
        }else if(best_pred(i,j) == S2_INC)
        {
            j--;
        }else
        {
            path.push_back({i-1,j-1});
            // if(i == 0 || j == 0) break;
            i--;
            j--;
        }
    }
    std::reverse(path.begin(), path.end());

    return path;
}

} // namespace mxm


#endif // _METRIC_STRING_H_
