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
    Vector<float> expected{0.1f, -0.6, 0.};

    if(!isZero(ret - expected))
    {
        std::cout << mxm::to_string(ret) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testBoundaryNormalVector01()
{
    Matrix<float> triangle({2,3},{1,0, 0,1, 0,0}, COL);
    Matrix<float> ret = splx::boundaryNormalVectors(triangle);
    Matrix<float> expected({2,3},{-1,0, 0,-1, mxm::sqrt(2.f)/2.f, mxm::sqrt(2.f)/2.f}, COL);
    if(!isZero(ret - expected, nullptr, 1e-5))
    {
        std::cout << mxm::to_string(ret) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testBoundaryNormalVector02()
{
    Matrix<float> tetra({3,4},{1,0,0, 0,1,0, 0,0,1, 0,0,0}, COL);
    Matrix<float> ret = splx::boundaryNormalVectors(tetra);
    Matrix<float> expected({3,4},{-1,0,0, 0,-1,0, 0,0,-1,  mxm::sqrt(3.f)/3.f, mxm::sqrt(3.f)/3.f, mxm::sqrt(3.f)/3.f}, COL);
    if(!isZero(ret - expected, nullptr, 1e-5))
    {
        std::cout << mxm::to_string(ret) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testSimplexDistanceToPoint()
{
    Matrix<float> triangle({3,3},{1,0,0, 0,1,0, 0,0,1}, COL);
    Matrix<float> pts( fixRow(3) ,{0.f,0,0, 1,1,1, 1,0,0}, COL);
    Vector<float> ref_pt{0.f,0,0};
    // auto dist = splx::distanceToPoints(triangle, pts, &ref_pt);
    auto ret = splx::distanceBoundaryToPoints(triangle, pts);
    Matrix<float> expected(fixRow(3), {
        -mxm::sqrt(6.f)/6, -mxm::sqrt(6.f)/6, -mxm::sqrt(6.f)/6,
        -mxm::sqrt(6.f)/6, -mxm::sqrt(6.f)/6, -mxm::sqrt(6.f)/6,
        -mxm::sqrt(6.f)/2, 0, 0}, COL);
    if(!isZero(ret - expected, nullptr, 1e-5))
    {
        std::cout << mxm::to_string(ret) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testSimplexDistanceToPoint02()
{
    size_t dim = 3;
    // Matrix<float> simplex(fixRow(dim),{1,0,0, 0,1,0}, COL);
    Matrix<float> simplex(fixRow(dim),{1,0,0, 0,1,0}, COL);
    Matrix<float> pts( fixRow(dim) ,{0.f,0,0}, COL);
    Vector<float> ref_pt{0.f,0,1};
    auto dist = splx::distanceSubspaceToPoints(simplex, pts);
    // auto dist_bdry = splx::distanceBoundaryToPoints(simplex, pts);
    Vector<float> expected{sqrt(0.5f)};
    if(!isZero(dist - expected, nullptr, 1e-5))
    {
        std::cout << mxm::to_string(dist) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testSimplex()
{
    testAffineCoordinate01();
    testBarycentricCoordinate01();
    testLebesgueMeasure();
    testPointsInsideSimplex();
    testBoundaryNormalVector01();
    testBoundaryNormalVector02();
    testSimplexDistanceToPoint();
    testSimplexDistanceToPoint02();
}

} // namespace mxm

#endif // _TEST_SIMPLEX_H_
