#if !defined(_TEST_BVH_H)
#define _TEST_BVH_H

#include "test_config.h"

#if TEST_AVAILABLE_ALL

#include "mxm/spatial_bvh.h"
#include "mxm/random.h"
#include <memory>
#include <set>

namespace mxm{

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
        2.000000, 0.000000, 0.000000, 6.000000, 1.000000, 0.000000}, COL);

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
    0.5, 0.5, -0.5}, COL));

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
        23, 22, 21}, COL));

    bvh::PrimitiveMeshTree tree(vertex_buffer, index_buffer);
    tree.build(1, false);

    if(1){
        Ray<> ray({-3, 0.1, 0}, {1,0,0});
        auto records = tree.hit(ray, bvh::eMultiHit);
        if(records.size() != 2)
        {
            std::cout << "hit_cnt: " << records.size() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        for(const auto & record: records)
        {
            if( ! isZero(ray(record.t) - tree.primitive(record.prim_idx).matmul(record.coeff), nullptr, 1e-4))
            {
                std::cout
                    << "ray(record.t): " << mxm::to_string (ray(record.t))
                    << ", primitive: " << mxm::to_string (tree.primitive(record.prim_idx).matmul(record.coeff))
                    << std::endl;
                throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
            }
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
    std::shared_ptr<Mat> pts(new Mat(fixRow(2), {
        0.749941, 0.910648, 0.132996, 0.181847, 0.982361, 0.263803, 0.095355, 0.145539, 0.282673, 0.136069, 0.802111, 0.869292, 0.077557, 0.579705, 0.627384, 0.549860, 0.008094, 0.144955, 0.680287, 0.853031, 0.533933, 0.622055, 0.438667, 0.350952, 0.199551, 0.513250, 0.138001, 0.401808, 0.382333, 0.075967, 0.762421, 0.239916, 0.040471, 0.123319, 0.252956, 0.183908, 0.504771, 0.239953, 0.822605, 0.417267, 0.981723, 0.049654, 0.823455, 0.902716, 0.301827, 0.944787, 0.047944, 0.490864, 0.247848, 0.489253, 0.544056, 0.337719, 0.887726, 0.900054, 0.343840, 0.369247, 0.033269, 0.111203, 0.162872, 0.780252, 0.877364, 0.389739, 0.210302, 0.241691, 0.272753, 0.403912, 0.492442, 0.096455, 0.223420, 0.131973, 0.490301, 0.942051, 0.954943, 0.956135, 0.651969, 0.575209, 0.757504, 0.059780, 0.435246, 0.234780, 0.552881, 0.353159, 0.053153, 0.821194, 0.322511, 0.015403, 0.404981, 0.043024, 0.908434, 0.168990, 0.809137, 0.649116, 0.258296, 0.731722, 0.129847, 0.647746, 0.493327, 0.450924, 0.378500, 0.547009,
        0.718470, 0.296321, 0.287805, 0.744693, 0.623436, 0.188955, 0.811874, 0.686775, 0.312508, 0.183511, 0.385711, 0.368485, 0.343930, 0.625619, 0.815769, 0.780227, 0.669285, 0.081126, 0.379419, 0.929386, 0.458497, 0.775713, 0.303614, 0.486792, 0.910565, 0.435859, 0.488618, 0.446784, 0.755790, 0.306349, 0.108062, 0.508509, 0.393590, 0.510772, 0.870187, 0.817628, 0.360861, 0.794831, 0.917118, 0.644318, 0.543805, 0.378609, 0.140144, 0.811580, 0.199873, 0.532826, 0.948925, 0.350727, 0.990110, 0.939002, 0.240076, 0.875943, 0.016521, 0.550156, 0.388615, 0.622475, 0.779689, 0.587045, 0.476638, 0.207742, 0.540138, 0.301246, 0.018451, 0.470923, 0.900183, 0.230488, 0.194495, 0.844309, 0.885998, 0.194764, 0.441223, 0.225922, 0.147829, 0.170708, 0.239502, 0.227664, 0.799653, 0.435699, 0.473015, 0.311102, 0.089823, 0.923380, 0.644551, 0.430207, 0.633064, 0.184816, 0.584382, 0.904881, 0.726654, 0.979748, 0.354638, 0.438870, 0.680407, 0.111119, 0.707322, 0.258065, 0.162329, 0.408720, 0.133736, 0.594896
    }));

    bvh::PointCloudTree tree(pts);
    tree.build(4, false);

    Vec target_pt({.5, .5});
    size_t k = 5;

    auto result = tree.nearestNeighborSearch(target_pt, k);
    std::multimap<float, size_t> expected;

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
#else
void testBvh(){}
#endif
} //scope of namespace mxm
#endif // _TEST_BVH_H
