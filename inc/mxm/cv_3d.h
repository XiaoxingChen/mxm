#if !defined(_CV_3D_H)
#define _CV_3D_H

#include "cv_basic.h"
#include "rigid_transform.h"

namespace mxm
{

template<typename DType>
Rotation<DType>
icpFindRotation(
    const Matrix<DType>& centralized_pts1,
    const Matrix<DType>& centralized_pts2)
{
    auto mat_w = centralized_pts2.matmul(centralized_pts1.T());
    auto u_d_vt = svd(mat_w);
    Matrix<DType> mat_r = u_d_vt[0].matmul(u_d_vt[2]);

    auto rot = Rotation<DType>::fromMatrix(mat_r);
    return rot;
}



} // namespace mxm


#endif // _CV_3D_H
