#if !defined(_TEST_INTERP_H_)
#define _TEST_INTERP_H_

#include "test_config.h"

#if TEST_AVAILABLE_ALL

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

void testBilinearInterpolation()
{
    {
        Mat square({2,2},{1,2,2,3});
        Mat pos(fixRow(2), {0.5,0.5, 0.25,0.25, 0,0}, Mat::COL);

        Mat expected(fixCol(1), {2, 1.5, 1});

        auto ret = interp::bilinearUnitSquare(pos, square);
        if((ret-expected).norm() > eps())
        {
            std::cout << mxm::to_string(ret) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {

    }




}

void testInterpolation()
{
    testTriangular();
    testBilinearInterpolation();
}
#else
void testInterpolation(){}
#endif


#endif // _TEST_INTERP_H_
