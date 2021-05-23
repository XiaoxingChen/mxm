#if !defined(_TEST_RIGID_TRANSFORM_H)
#define _TEST_RIGID_TRANSFORM_H

#include "mxm/rigid_transform.h"
using namespace mxm;

inline void testRigidTransform()
{
    size_t dim = 3;
    RigidTrans pose(Vec::zeros(dim), Rotation::identity(dim));
    if((pose.asMatrix() - Mat::identity(dim + 1)).norm() > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}


#endif // _TEST_RIGID_TRANSFORM_H
