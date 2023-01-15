#if !defined(_GEOMETRY_SIMPLEX_H_)
#define _GEOMETRY_SIMPLEX_H_

#include "mxm/linalg_solve.h"

namespace mxm
{

namespace splx
{
#if 0
template<typename DType>
Matrix<DType> dimReduction(const Matrix<DType>& vertices)
{
    Matrix<DType> basis = vertices(Block({},{1, vertices.shape(1)}));
    basis -= vertices(Col(0));
    auto qr_result = qr::decomposeByRotation(basis);
    return qr_result[1];
}
#endif

// find the barycentric coordinates of points with respect to the simplex
// M: number of points
// N: Dimension
// pts_in: shape(N,M)
// simplex: shape(N,N)
// return: shape(N-1,M)
// If point is not in the N-1 subspace of simplex, 
// the projection of the point will be used for calculation instead. 
template<typename DType>
Matrix<DType> barycentricCoordinate(
    const Matrix<DType>& pts_in, 
    const Matrix<DType>& simplex)
{
    // dim reduction

    // move vertex_0 to center
    Matrix<DType> basis = simplex(Block({},{1, simplex.shape(1)}));
    basis -= simplex(Col(0));

    Matrix<DType> pts = pts_in - simplex(Col(0));

    // rotate
    auto qr_result = qr::decomposeByRotation(basis);
    basis = qr_result[1](Block({0, simplex.shape(0)-1}, {}));

    // std::cout << mxm::to_string(qr_result) << std::endl;

    // solve the coordinates
    pts = qr_result[0].T().matmul(pts);
    Matrix<DType> pts_proj = pts(Block({0, simplex.shape(0)-1}, {}));
    return qr::solve(basis, pts_proj);
}

} // namespace splx

} // namespace mxm

#endif // _GEOMETRY_SIMPLEX_H_
