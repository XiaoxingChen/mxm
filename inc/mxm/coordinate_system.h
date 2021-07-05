#if !defined(_COORDINATE_SYSTEM_H_)
#define _COORDINATE_SYSTEM_H_

#include "linalg.h"

namespace mxm
{

// input: x, y, z
// north pole: z axis
// return latitude, longitude, height(radius)
// Reference:
// [1] https://en.wikipedia.org/wiki/Spherical_coordinate_system
template<typename DType>
Matrix<DType> cartesianToSpherical(const Matrix<DType>& cartesian)
{
    Matrix<DType> spherical(cartesian.shape());
    for(size_t i = 0; i < cartesian.shape(1); i++)
    {
        const DType& x = cartesian(0, i);
        const DType& y = cartesian(1, i);
        const DType& z = cartesian(2, i);

        spherical(0, i) = atan2(z, sqrt(x*x + y*y));
        spherical(1, i) = atan2(y, x);
        spherical(2, i) = sqrt(x*x + y*y + z*z);
    }
    return spherical;
}

// input: latitude, longitude, height(radius)
// north pole: z axis
// return: x, y, z
// Reference:
// [1] https://en.wikipedia.org/wiki/Spherical_coordinate_system
template<typename DType>
Matrix<DType> sphericalToCartesian(const Matrix<DType>& spherical)
{
    Matrix<DType> cartesian(spherical.shape());
    for(size_t i = 0; i < cartesian.shape(1); i++)
    {
        const DType& lon = spherical(0, i);
        const DType& lat = spherical(1, i);
        const DType& r   = spherical(2, i);

        cartesian(0, i) = r * cos(lon) * cos(lat);
        cartesian(1, i) = r * cos(lon) * sin(lat);
        cartesian(2, i) = r * sin(lon);
    }
    return cartesian;
}

} // namespace mxm


#endif // _COORDINATE_SYSTEM_H_
