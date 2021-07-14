#if !defined(_TEST_RIGID_TRANSFORM_H)
#define _TEST_RIGID_TRANSFORM_H
#include "test_config.h"

#if TEST_AVAILABLE_RIGID_TRANSFORM
#include "mxm/rigid_transform.h"
#include "mxm/rotation.h"
using namespace mxm;

inline void testRigidTransform()
{
    size_t dim = 3;
    RigidTrans pose(Vec::zeros(dim), Rotation<float>::identity(dim));
    if((pose.asMatrix() - Mat::identity(dim + 1)).norm() > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}
#else
inline void testRigidTransform(){}
#endif


#endif // _TEST_RIGID_TRANSFORM_H
