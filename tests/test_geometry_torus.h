#if !defined(__TEST_GEOMETRY_TORUS_H__)
#define __TEST_GEOMETRY_TORUS_H__

#include "mxm/geometry_torus.h"

using namespace mxm;

inline void testGenerateCircle()
{
    Matrix<float> vertex_buffer;
    Matrix<size_t> index_buffer;

    vertex_buffer = generateNTorus<2, float>({1.f}, {4});
    if(isZero(vertex_buffer - Matrix<float>({2,4},{1,0, 0,1, -1,0, 0,-1}, COL)))
    {
        std::cout << mxm::to_string(vertex_buffer) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    vertex_buffer = generateNTorus<3, float>({0.f, 1}, {1, 4});
    if(isZero(vertex_buffer - Matrix<float>({3,4},{0,1,0, 0,0,1, 0,-1,0, 0,0,-1}, COL)))
    {
        std::cout << mxm::to_string(vertex_buffer) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testGeometryTorus()
{
    testGenerateCircle();
}

#endif // __TEST_GEOMETRY_TORUS_H__
