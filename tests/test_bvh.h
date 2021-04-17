#if !defined(_TEST_BVH_H)
#define _TEST_BVH_H

#include "mxm/spatial_bvh.h"
#include <memory>

using namespace mxm;

inline std::shared_ptr<Mat> createTriangleBand(
    size_t half_n,
    Matrix<size_t>& indices)
{
    size_t dim = 3;
    std::shared_ptr<Mat> ret(new Mat({dim, half_n * 2}));
    Mat& vertices = *ret;
    // indices.clear();
    indices = Matrix<size_t>({3, half_n});
    // std::vector<std::vector< size_t>> indices;
    for(size_t n = 0; n < half_n; n++)
    {
        size_t i = 2*n;
        vertices(Col(i)) = Vec({static_cast<FloatType>(i), 0, 0});
        vertices(Col(i + 1)) = Vec({static_cast<FloatType>(i), 1, 0});

        if(i < half_n - 1)
        {
            indices(Col(i)) = Vector<size_t>({i, i + 3, i + 1});
            indices(Col(i+1)) = Vector<size_t>({i, i + 2, i + 3});
            // indices.push_back({i, i + 3, i + 1});
            // indices.push_back({i, i + 2, i + 3});
        }
    }
    return ret;
}

inline void testSort()
{

}

inline void testBuildTree()
{


}

inline void testBuildTree2()
{
    std::cout << "testBuildTree2()" << std::endl;
    size_t dim = 3;
    std::shared_ptr<Matrix<size_t>> index_buffer(new Matrix<size_t>);

    size_t triangle_num = 40;
    auto vertex_buffer = createTriangleBand(40, *index_buffer);

    bvh::PrimitiveMeshTree tree(vertex_buffer, index_buffer);
    tree.build(4, false);
}

inline void testMultiHit()
{
    std::shared_ptr<Mat> vertex_buffer(new Mat({3, 24},{
    -0.5, -0.5, 0.5,
    0.5, -0.5, 0.5,
    -0.5, 0.5, 0.5,
    0.5, 0.5, 0.5,
    0.5, -0.5, 0.5,
    -0.5, -0.5, 0.5,
    0.5, -0.5, -0.5,
    -0.5, -0.5, -0.5,
    0.5, 0.5, 0.5,
    0.5, -0.5, 0.5,
    0.5, 0.5, -0.5,
    0.5, -0.5, -0.5,
    -0.5, 0.5, 0.5,
    0.5, 0.5, 0.5,
    -0.5, 0.5, -0.5,
    0.5, 0.5, -0.5,
    -0.5, -0.5, 0.5,
    -0.5, 0.5, 0.5,
    -0.5, -0.5, -0.5,
    -0.5, 0.5, -0.5,
    -0.5, -0.5, -0.5,
    -0.5, 0.5, -0.5,
    0.5, -0.5, -0.5,
    0.5, 0.5, -0.5}, Mat::COL));

    std::shared_ptr<Matrix<size_t>> index_buffer( new Matrix<size_t>({3,12},{
        0 , 1 , 2 ,
        3 , 2 , 1 ,
        4 , 5 , 6 ,
        7 , 6 , 5 ,
        8 , 9 , 10,
        11, 10, 9 ,
        12, 13, 14,
        15, 14, 13,
        16, 17, 18,
        19, 18, 17,
        20, 21, 22,
        23, 22, 21}, Mat::COL));

    bvh::PrimitiveMeshTree tree(vertex_buffer, index_buffer);
    tree.build(1, false);

    if(1){
        Ray ray({-3, 0.1, 0}, {1,0,0});
        auto hit_cnt = tree.multiHit(ray);
        if(hit_cnt != 2)
        {
            std::cout << "hit_cnt: " << hit_cnt << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }


}

void testBvh()
{
    testSort();
    testBuildTree();
    testBuildTree2();
    testMultiHit();
}

#endif // _TEST_BVH_H
