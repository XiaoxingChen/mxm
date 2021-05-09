#if !defined(_TEST_BVH_H)
#define _TEST_BVH_H

#include "mxm/spatial_bvh.h"
#include "mxm/random.h"
#include <memory>
#include <set>

using namespace mxm;

// A triangle band from [0,half_n]
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

    const auto & node_buffer = tree.nodeBuffer();

    Mat expected({6, 31},
        {0.000000, 0.000000, 0.000000, 40.000000, 1.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 20.000000, 1.000000, 0.000000,
        20.000000, 0.000000, 0.000000, 40.000000, 1.000000, 0.000000,
        20.000000, 0.000000, 0.000000, 30.000000, 1.000000, 0.000000,
        30.000000, 0.000000, 0.000000, 40.000000, 1.000000, 0.000000,
        30.000000, 0.000000, 0.000000, 36.000000, 1.000000, 0.000000,
        34.000000, 0.000000, 0.000000, 40.000000, 1.000000, 0.000000,
        34.000000, 0.000000, 0.000000, 38.000000, 1.000000, 0.000000,
        36.000000, 0.000000, 0.000000, 40.000000, 1.000000, 0.000000,
        30.000000, 0.000000, 0.000000, 32.000000, 1.000000, 0.000000,
        32.000000, 0.000000, 0.000000, 36.000000, 1.000000, 0.000000,
        20.000000, 0.000000, 0.000000, 26.000000, 1.000000, 0.000000,
        24.000000, 0.000000, 0.000000, 30.000000, 1.000000, 0.000000,
        24.000000, 0.000000, 0.000000, 28.000000, 1.000000, 0.000000,
        26.000000, 0.000000, 0.000000, 30.000000, 1.000000, 0.000000,
        20.000000, 0.000000, 0.000000, 22.000000, 1.000000, 0.000000,
        22.000000, 0.000000, 0.000000, 26.000000, 1.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 10.000000, 1.000000, 0.000000,
        10.000000, 0.000000, 0.000000, 20.000000, 1.000000, 0.000000,
        10.000000, 0.000000, 0.000000, 16.000000, 1.000000, 0.000000,
        14.000000, 0.000000, 0.000000, 20.000000, 1.000000, 0.000000,
        14.000000, 0.000000, 0.000000, 18.000000, 1.000000, 0.000000,
        16.000000, 0.000000, 0.000000, 20.000000, 1.000000, 0.000000,
        10.000000, 0.000000, 0.000000, 12.000000, 1.000000, 0.000000,
        12.000000, 0.000000, 0.000000, 16.000000, 1.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 6.000000, 1.000000, 0.000000,
        4.000000, 0.000000, 0.000000, 10.000000, 1.000000, 0.000000,
        4.000000, 0.000000, 0.000000, 8.000000, 1.000000, 0.000000,
        6.000000, 0.000000, 0.000000, 10.000000, 1.000000, 0.000000,
        0.000000, 0.000000, 0.000000, 2.000000, 1.000000, 0.000000,
        2.000000, 0.000000, 0.000000, 6.000000, 1.000000, 0.000000}, Mat::COL);

    for(size_t i = 0; i < node_buffer.size(); i++)
    {
        if((node_buffer.at(i).aabb.min() - expected(Block({0,3},{i, i+1}))).norm() > eps())
        {
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
        if((node_buffer.at(i).aabb.max() - expected(Block({3,6},{i, i+1}))).norm() > eps())
        {
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

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

void testRadiusSearch()
{
    size_t dim = 2;
    std::shared_ptr<Mat> pts(new Mat(random::uniform<FloatType>({dim, 100})));

    bvh::PointCloudTree tree(pts);
    tree.build(4, false);

    Vec target_pt({.5, .5});
    FloatType radius = 0.5;

    auto result = tree.radiusSearch(target_pt, radius);
    std::multimap<FloatType, size_t> expected;
    for(size_t idx = 0; idx < pts->shape(1); idx++)
    {
        auto distance = (target_pt - (*pts)(Col(idx))).norm();
        if(distance < radius)
            expected.insert({distance , idx});
    }

    if(expected != result)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

void testNearestNeighborSearch()
{
    size_t dim = 2;
    std::shared_ptr<Mat> pts(new Mat(random::uniform<FloatType>({dim, 100})));

    bvh::PointCloudTree tree(pts);
    tree.build(4, false);

    Vec target_pt({.5, .5});
    size_t k = 5;

    auto result = tree.nearestNeighborSearch(target_pt, k);
    std::multimap<FloatType, size_t> expected;

    for(size_t i = 0; i < pts->shape(1); i++)
    {
        expected.insert({((*pts)(Col(i)) - target_pt).norm(), i});
    }
    while(expected.size() > k) expected.erase(std::prev(expected.end()));

    if(expected != result)
    {
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

void testBvh()
{
    testSort();
    testBuildTree();
    testBuildTree2();
    testMultiHit();
    testRadiusSearch();
    testNearestNeighborSearch();
}

#endif // _TEST_BVH_H
