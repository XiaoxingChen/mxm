#if !defined(_TEST_INTERP_H_)
#define _TEST_INTERP_H_

#include <iostream>
#include "mxm/interpolation.h"

void testTriangular()
{
    Mat tex_coord({2,3}, {0,0, 1,0, 0.5, 0.5*sqrt(3.)}, Mat::COL);
    auto ret = interp::triangular(Vec({.5, .5}), tex_coord, tex_coord);
    if((ret - Vec::ones(2) * .5).norm() > eps())
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

void testInterpolation()
{
    testTriangular();
}



#endif // _TEST_INTERP_H_
