#if !defined(_TEST_SIMPLEX_H_)
#define _TEST_SIMPLEX_H_

#include "mxm/geometry_simplex.h"

namespace mxm
{

inline void testBarycentricCoordinate01()
{
    Matrix<float> pts(fixRow(3), {
        1,0,0, 
        0,1,0, 
        0,0.5,0.5}, COL);
    Matrix<float> triangle(fixRow(3), {1,0,0, 0,1,0, 0,0,1}, COL);
    Matrix<float> coord = splx::barycentricCoordinate(pts, triangle);

    Matrix<float> expected(fixRow(2), {
        0,0, 
        1,0, 
        0.5,0.5}, COL);
    if(!isZero(coord - expected, nullptr, 1e-7))
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testSimplex()
{
    testBarycentricCoordinate01();
}

} // namespace mxm

#endif // _TEST_SIMPLEX_H_
