#if !defined(__TEST_GRID_MAP_H__)
#define __TEST_GRID_MAP_H__

#include "test_config.h"

#if TEST_AVAILABLE_ALL
#include "mxm/spatial_grid_map.h"
#endif

namespace mxm{

#if TEST_AVAILABLE_ALL
void testQuantize01()
{
    Mat pts({2,3}, {3,3, 9,9, 20,4}, COL);
    Matrix<uint32_t> quantized_pts = quantize(pts, 128);
    Matrix<uint32_t> expected({2,3}, {0,45,127, 0,45,7});
    if(quantized_pts != expected)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testGridMap()
{
    testQuantize01();
}
#else
void testGridMap(){}
#endif

} //scope of namespace mxm
#endif // __TEST_GRID_MAP_H__
