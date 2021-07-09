#if !defined(_TEST_COORDINATE_SYSTEM_H_)
#define _TEST_COORDINATE_SYSTEM_H_

#include "mxm/coordinate_system.h"

namespace mxm
{
inline void testCoordinateSystem()
{
    {
        Matrix<float> cartesian(fixRow(3),{
        1,0,0,
        1,1,1,
        0,0,1
        }, COL);

        auto inv_inv = cartesianFromSpherical(sphericalFromCartesian(cartesian));

        if(norm(inv_inv - cartesian) > std::numeric_limits<float>::epsilon())
        {
            std::cout << mxm::to_string(inv_inv) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }
}

} // namespace mxm



#endif // _TEST_COORDINATE_SYSTEM_H_
