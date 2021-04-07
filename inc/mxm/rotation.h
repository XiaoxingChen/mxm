#if !defined(_ROTATION_H_)
#define _ROTATION_H_

#include "low_dimensional_rotation.h"
#include "full_dimensional_rotation.h"



namespace mxm
{


using Rotation = FullDimensionalRotation<FloatType>;
// using Rotation = LowDimensionalRotation;
} // namespace mxm


#endif // _ROTATION_H_
