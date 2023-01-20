#if !defined(_TEST_SIMPLEX_H_)
#define _TEST_SIMPLEX_H_

#include "mxm/geometry_simplex.h"
#include "mxm/interpolation.h"

namespace mxm
{

inline void testAffineCoordinate01()
{
    Matrix<float> pts(fixRow(3), {
        1,0,0, 
        0,1,0, 
        0,0.5,0.5}, COL);
    Matrix<float> triangle(fixRow(3), {1,0,0, 0,1,0, 0,0,1}, COL);
    Matrix<float> coord = splx::affineCoordinate(pts, triangle);

    Matrix<float> expected(fixRow(2), {
        0,0, 
        1,0, 
        0.5,0.5}, COL);
    if(!isZero(coord - expected, nullptr, 1e-7))
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testBarycentricCoordinate01()
{
    Mat tex_coord({2,3}, {0,0, 1,0, 0.5, 0.5f*sqrt(3.f)}, COL);
    auto pts = Vec({.5, .5});
    auto ret1 = interp::triangular(pts, tex_coord, tex_coord);

    Matrix<float> weights = splx::barycentricCoordinate(pts, tex_coord);
    auto ret2 = tex_coord.matmul(weights);
    // std::cout << mxm::to_string(ret1) << std::endl;
    // std::cout << mxm::to_string(ret2) << std::endl;

    if(!isZero(ret1 - ret2, nullptr, 1e-7))
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testLebesgueMeasure()
{
    Matrix<float> triangle({3,3},{1,0,0, 0,1,0, 0,0,0}, COL);
    float area = splx::lebesgueMeasure(triangle);
    if(fabs(area - 0.5) > eps<float>())
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    Matrix<float> tetrahedron({3,4},{0,0,0, 1,0,0, 0,1,0, 0,0,1}, COL);
    float volume = splx::lebesgueMeasure(tetrahedron);
    if(fabs(volume - 1./6) > eps<float>())
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testPointsInsideSimplex()
{
    Matrix<float> triangle({2,3},{1,0, 0,1, 0,0}, COL);
    Matrix<float> pts(fixRow(2), {0.1,0.1, 0.5,1.1, 0.5,0.5}, COL);
    auto ret = splx::arePointsInside(pts, triangle);
    Vector<int8_t> expected{1,0,1};
    
    if(!isZero(ret - expected))
    {
        std::cout << mxm::to_string(ret) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testSimplex()
{
    testAffineCoordinate01();
    testBarycentricCoordinate01();
    testLebesgueMeasure();
    testPointsInsideSimplex();
}

} // namespace mxm

#endif // _TEST_SIMPLEX_H_
