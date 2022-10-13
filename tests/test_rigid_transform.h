#if !defined(_TEST_RIGID_TRANSFORM_H)
#define _TEST_RIGID_TRANSFORM_H
#include "test_config.h"

#if TEST_AVAILABLE_RIGID_TRANSFORM
#include "mxm/rigid_transform.h"
#include "mxm/rotation.h"
namespace mxm{

inline void testRigidTransform()
{
    const size_t DIM = 3;
    RigidTrans pose(Vec::zeros(DIM), Rotation<float, DIM>::identity());
    if((pose.asMatrix() - Mat::identity(DIM + 1)).norm() > eps())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
}
#else
inline void testRigidTransform(){}
#endif

} //scope of namespace mxm
#endif // _TEST_RIGID_TRANSFORM_H
