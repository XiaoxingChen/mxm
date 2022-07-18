#if !defined(__TEST_MODEL_LIDAR_H__)
#define __TEST_MODEL_LIDAR_H__

#include "mxm/model_lidar.h"

using namespace mxm;

inline void testModelLidar()
{
    LidarRotating<float, 2> lidar_2d;
    LidarRotating<float, 3> lidar_3d;

    lidar_2d.setHorizontalResolution(4);

    lidar_3d.setVerticalAngles({-1, 0, 1}).setHorizontalResolution(4);

    auto dirs_2d = lidar_2d.castRayDirection();
    auto dirs_3d = lidar_3d.castRayDirection();

    // std::cout << mxm::to_string(dirs_2d) << std::endl;
    // std::cout << mxm::to_string(dirs_3d) << std::endl;
}

#endif // __TEST_MODEL_LIDAR_H__
