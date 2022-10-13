#if !defined(__TEST_BSP_TREE_H__)
#define __TEST_BSP_TREE_H__

#include "test_config.h"

#if TEST_AVAILABLE_ALL
#include "mxm/spatial_bsp.h"

namespace mxm{

inline void testAxisPartition()
{
    {
        Mat points({2, 4}, {0,0, 0,1, 1,0, 1,1}, COL);
        std::vector<size_t> index_buffer{3,2,1,0};
        Vec thresholds({0.5, 0,5});

        auto p = bsp::axisPartition(index_buffer, {0,4}, points, 1, thresholds);

        if(p != 2)
        {
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
        if(Vector<size_t>(index_buffer) != Vector<size_t>({0,2,1,3}))
        {
            std::cout << mxm::to_string(index_buffer) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {
        Mat points({2, 4}, {0,0, 0,1, 1,0, 1,1}, COL);
        std::vector<size_t> index_buffer{3,2,1,0};
        Vec thresholds({0.5, 0,5});

        auto partitioners = bsp::fullAxesPartition(index_buffer, {0,4}, points, thresholds);


        if(Vector<size_t>(index_buffer) != Vector<size_t>({0,1,2,3}))
        {
            std::cout << mxm::to_string(index_buffer) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
        if(Vector<size_t>(partitioners) != Vector<size_t>({0,1,2,3,4}))
        {
            std::cout << mxm::to_string(partitioners) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    if(0){
        std::shared_ptr<Mat> p_points(new Mat ({2, 4}, {0,0, 0,1, 1,0, 1,1}, COL));
        bsp::PointCloudTree tree(p_points);
        tree.build(1, true);
        for(const auto & node: tree.nodeBuffer())
        {
            std::cout << node.width << std::endl;
        }
    }

    if(0){
        std::shared_ptr<Mat> p_points(new Mat ({3, 9}));
        for(size_t i = 0; i < 8; i++) (*p_points)(Col(i)) = binaryToVector<FloatType>(3, i);
        (*p_points)(Col(8)) = Vec({99,99,99});

        bsp::PointCloudTree tree(p_points);
        tree.build(1, true);
        for(const auto & node: tree.nodeBuffer())
        {
            std::cout << node.width << std::endl;
        }
    }
}


inline void testBspTree()
{
    testAxisPartition();
}
#else
inline void testBspTree(){}
#endif

} //scope of namespace mxm
#endif // __TEST_BSP_TREE_H__
