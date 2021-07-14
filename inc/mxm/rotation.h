#if !defined(_ROTATION_H_)
#define _ROTATION_H_

#include "low_dimensional_rotation.h"
#include "full_dimensional_rotation.h"



namespace mxm
{


template<typename DType>
using Rotation = FullDimensionalRotation<DType>;
// using Rotation = LowDimensionalRotation;
} // namespace mxm


#endif // _ROTATION_H_
