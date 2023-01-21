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

// find the affine coordinates of points with respect to the simplex
// K: number of points
// M: Cartisan Dimension
// N: Simplex Dimension
// pts_in: shape(M,K)
// simplex: shape(M,N+1)
// return: shape(N,K)
// If point is not in the N subspace of simplex, 
// the projection of the point will be used for calculation instead. 
// reference: 
// [1] https://en.wikipedia.org/wiki/Affine_space#Coordinates
template<typename DType>
Matrix<DType> affineCoordinate(
    const Matrix<DType>& pts_in, 
    const Matrix<DType>& simplex)
{
    // dim reduction
    size_t simplex_dim = simplex.shape(1) - 1;

    // move vertex_0 to center
    Matrix<DType> basis = simplex(Block({},{1, simplex.shape(1)}));
    basis -= simplex(Col(0));

    Matrix<DType> pts = pts_in - simplex(Col(0));

    // rotate
    auto qr_result = qr::decomposeByRotation(basis);
    basis = qr_result[1](Block({0, simplex_dim}, {}));

    // std::cout << mxm::to_string(qr_result) << std::endl;

    // solve the coordinates
    pts = qr_result[0].T().matmul(pts);
    Matrix<DType> pts_proj = pts(Block({0, simplex_dim}, {}));
    return qr::solve(basis, pts_proj);
}

// find the barycentric coordinates of points with respect to the simplex
// K: number of points
// M: Cartisan Dimension
// N: Simplex Dimension
// pts_in: shape(M,K)
// simplex: shape(M,N+1)
// return: shape(N+1,K)
template<typename DType>
Matrix<DType> barycentricCoordinate(
    const Matrix<DType>& pts_in, 
    const Matrix<DType>& simplex)
{
    Matrix<DType> affine = affineCoordinate(pts_in, simplex);
    // Matrix<DType> bary({affine.shape(0) + 1, affine.shape(1)});
    Matrix<DType> row0 = mxm::sum(affine, 0) * DType(-1.) + DType(1.);
    return vstack(row0, affine);
}

template<typename DType>
DType lebesgueMeasure(const Matrix<DType>& simplex)
{
    Matrix<DType> basis = simplex(Block({},{1, simplex.shape(1)}));
    basis -= simplex(Col(0));
    auto qr_result = qr::decomposeByRotation(basis);
    size_t simplex_dim = simplex.shape(1) - 1;
    DType ret(1.);
    for(size_t i = 0; i < simplex_dim; i++)
    {
        ret *= qr_result[1](i,i);
    }
    DType alpha(1);
    for(size_t i = 2; i <= simplex_dim; i++) alpha /= DType(i);
    ret *= alpha;

    return ret;
}

template<typename DType>
Vector<int8_t> arePointsInside(
    const Matrix<DType>& pts_in, 
    const Matrix<DType>& simplex,
    DType tol=DType(0.))
{
    assert(simplex.shape(0) + 1 == simplex.shape(1) && "simplex shape mismatch");
    assert(simplex.shape(0) == pts_in.shape(0) && "simplex shape mismatch");
    
    auto bary_coord = barycentricCoordinate(pts_in, simplex);
    Vector<int8_t> ret = Vector<int8_t>::ones(pts_in.shape(1));
    for(size_t j = 0; j < bary_coord.shape(1); j++)
    {
        for(size_t i = 0; i < bary_coord.shape(0); i++)
        {
            if(bary_coord(i,j) < DType(0.) - tol || bary_coord(i,j) > DType(1.) + tol)
            {
                ret(j) = 0;
                continue;
            }
        }
    }
    return ret;
}

template<typename DType>
Vector<DType> centroid(const Matrix<DType>& simplex)
{
    return sum(simplex, 1) / DType(simplex.shape(1)); 
}

#if 0
// reference:
// https://math.stackexchange.com/questions/4056099/circumcenter-of-the-n-simplex
template<typename DType>
Vector<DType> circumCenter(const Matrix<DType>& simplex, DType* radius=nullptr)
{
    
}
#endif

} // namespace splx

} // namespace mxm

#endif // _GEOMETRY_SIMPLEX_H_
