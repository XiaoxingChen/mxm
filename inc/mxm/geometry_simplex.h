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

// parameters
// N points, DIM dimension, M simplex dimension
// pts_in: shape(DIM, N)
// simplex: shape(DIM, M)
// return: shape(N, 1).
//   ret[i] > 0: inside
//   ret[i] == 0: on the boundary
//   ret[i] < 0: outside
template<typename DType>
Vector<DType> arePointsInside(
    const Matrix<DType>& pts_in,
    const Matrix<DType>& simplex)
{
    assert(simplex.shape(0) + 1 == simplex.shape(1) && "simplex shape mismatch");
    assert(simplex.shape(0) == pts_in.shape(0) && "simplex shape mismatch");

    auto bary_coord = barycentricCoordinate(pts_in, simplex);
    Vector<DType> ret(pts_in.shape(1));
    for(size_t j = 0; j < bary_coord.shape(1); j++)
    {
        ret(j) = mxm::min(bary_coord(Col(j)));
    }
    return ret;
}

template<typename DType>
Matrix<DType> boundaryNormalVectors(const Matrix<DType>& simplex)
{
    Matrix<DType> affine_basis({simplex.shape(0), simplex.shape(1) - 1});
    Matrix<DType> boundary_space({simplex.shape(0), simplex.shape(1) - 2});
    Matrix<DType> affine_boundary({simplex.shape(1) - 2, simplex.shape(1) - 1});
    Matrix<DType> ret({simplex.shape(0), simplex.shape(1)});
    for(size_t src_vtx_idx = 0; src_vtx_idx < simplex.shape(1); src_vtx_idx++)
    {
        size_t basis_idx = 0;
        std::vector<size_t> valid_dst;
        for(size_t dst_vtx_idx = 0; dst_vtx_idx < simplex.shape(1); dst_vtx_idx++)
        {
            if(dst_vtx_idx == src_vtx_idx) continue;
            affine_basis(Col(basis_idx)) = simplex(Col(dst_vtx_idx)) - simplex(Col(src_vtx_idx));
            valid_dst.push_back(dst_vtx_idx);
            basis_idx++;
        }
        for(size_t i = 0; i < valid_dst.size() - 1; i++)
        {
            boundary_space(Col(i)) = simplex(Col(valid_dst.at(i + 1))) - simplex(Col(valid_dst.at(i)));
        }

        affine_boundary = boundary_space.T().matmul(affine_basis);
        Matrix<DType> weights = orthogonalComplement(affine_boundary);
        if(sum(weights) < DType(0)) weights *= DType(-1);

        Matrix<DType> norm_v = affine_basis.matmul(weights.T());
        norm_v /= norm(norm_v);

        ret(Col(src_vtx_idx)) = norm_v;
    }
    return ret;
}

template<typename DType>
Vector<int8_t> normalVectorLeftUpRulePolarity(const Matrix<DType>& vectors)
{
    Vector<int8_t> ret = Vector<int8_t>::zeros(vectors.shape(1));
    for(size_t v_idx = 0; v_idx < vectors.shape(1); v_idx++)
    {
        for(size_t dim_idx = 0; dim_idx < vectors.shape(0); dim_idx++)
        {
            if(vectors(dim_idx, v_idx) > 0)
            {
                ret(v_idx) = 1;
                break;
            }
            if(vectors(dim_idx, v_idx) < 0)
            {
                ret(v_idx) = -1;
                break;
            }
        }
    }
    return ret;
}

#if 0
template<typename DType>
Vector<int8_t> arePixelsInside(
    const Matrix<DType>& bary_coord,
    const Matrix<DType>& simplex,
    DType tol)
{
    assert(simplex.shape(0) + 1 == simplex.shape(1) && "simplex shape mismatch");
    assert(simplex.shape(0) == bary_coord.shape(0) && "simplex shape mismatch");

    Vector<int8_t> ret = Vector<int8_t>::ones(bary_coord.shape(1));
    auto n_vec = boundaryNormalVectors(simplex);
    auto polarity = normalVectorLeftUpRulePolarity(n_vec);

    std::vector<size_t> exceed_idx;
    for(size_t j = 0; j < bary_coord.shape(1); j++)
    {
        exceed_idx.clear();
        for(size_t i = 0; i < bary_coord.shape(0); i++)
        {

            if(bary_coord(i,j) < -tol)
            {
                ret(j) = 0;
                continue;
            }
            if(bary_coord(i,j) < tol && polarity(i) < 0)
            {
                ret(j) = 0;
                continue;
            }
        }
    }
    return ret;
}
#endif

// parameters
// N points, DIM dimension, M simplex dimension
// pts_in: shape(DIM, N)
// simplex: shape(DIM, M)
// p_ref: reference point that specifies the negative side of the space
// return: shape(N)
template<typename DType>
Vector<DType>
distanceSubspaceToPoints(const Matrix<DType>& simplex, const Matrix<DType>& pts_in)
{
    // dim reduction
    size_t simplex_dim = simplex.shape(1) - 1;

    // move vertex_0 to center
    Matrix<DType> affine_basis = simplex(Block({},{1, end()})) - simplex(Col(0));

    Matrix<DType> pts = pts_in - simplex(Col(0));

    // rotate
    auto qr_result = qr::decomposeByRotation(affine_basis);

    // pts under rotated coordinate
    pts = qr_result[0].T().matmul(pts);


    Vector<DType> ret(pts_in.shape(1));

    for(size_t pt_idx = 0; pt_idx < pts_in.shape(1); pt_idx++)
    {
        ret(pt_idx) = norm(pts(Block({simplex_dim, end()}, {pt_idx, pt_idx+1})));
    }
    return ret;
}

// given a line as a simplex (DIM, 2)
// find the distance between the projected points and the line
// projected points are the points that projected to the subspace spanned from simplex and the origin point
template<typename DType>
Vector<DType>
distanceToPointsWithinSubspace(const Matrix<DType>& simplex, const Matrix<DType>& pts_in, const Vector<DType>& ref_pt)
{
    Matrix<DType> affine_basis = hstack(
        simplex(Block({},{1,end()})) - simplex(Col(0)),
        ref_pt - simplex(Col(0))
        );

    Matrix<DType> pts = pts_in - simplex(Col(0));

    // rotate
    auto qr_result = qr::decomposeByRotation(affine_basis);

    // pts under rotated coordinate
    pts = qr_result[0].T().matmul(pts);
    // auto projected_origin = qr_result[0].T().matmul(ref_pt);

    size_t dist_axis = simplex.shape(1) - 1;

    Vector<DType> ret = pts(Row(dist_axis));
    if(qr_result[1](dist_axis, dist_axis) > 0) ret *= DType(-1);
    return ret;
}

template<typename DType>
Matrix<DType>
distanceBoundaryToPoints(const Matrix<DType>& simplex, const Matrix<DType>& pts_in)
{
    Matrix<DType> ret({simplex.shape(1), pts_in.shape(1)});
    Matrix<DType> boundary_simplex({simplex.shape(0), simplex.shape(1) - 1});
    for(size_t vtx_idx = 0; vtx_idx < simplex.shape(1); vtx_idx++)
    {
        size_t boundary_vtx_idx = 0;
        for(size_t i = 0; i < simplex.shape(1); i++)
        {
            if(i == vtx_idx) continue;
            boundary_simplex(Col(boundary_vtx_idx)) = simplex(Col(i));
            boundary_vtx_idx++;
        }
        Vector<DType> ref_point = simplex(Col(vtx_idx));
        ret(Row(vtx_idx)) = distanceToPointsWithinSubspace(boundary_simplex, pts_in, ref_point).T();
    }
    return ret;
}

template<typename DType>
Vector<int8_t> arePointsInside(
    const Matrix<DType>& pts_in,
    const Matrix<DType>& simplex,
    DType boundary_width)
{
    assert(simplex.shape(0) + 1 == simplex.shape(1) && "simplex shape mismatch");
    assert(simplex.shape(0) == pts_in.shape(0) && "simplex shape mismatch");

    Vector<int8_t> ret = Vector<int8_t>::ones(pts_in.shape(1));
    auto n_vec = boundaryNormalVectors(simplex);
    auto polarity = normalVectorLeftUpRulePolarity(n_vec);
    auto distance_to_boundary = distanceBoundaryToPoints(simplex, pts_in);

    for(size_t j = 0; j < distance_to_boundary.shape(1); j++)
    {
        for(size_t vtx_idx = 0; vtx_idx < distance_to_boundary.shape(0); vtx_idx++)
        {
            if(distance_to_boundary(vtx_idx,j) > boundary_width)
            {
                ret(j) = 0;
                continue;
            }
            if(distance_to_boundary(vtx_idx,j) > -boundary_width && polarity(vtx_idx) < 0)
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

// reference:
// docs/simplex_geometry.html#h1_intersection_under_affine_frame
template<typename DType>
Vector<DType> intersect(
    const Matrix<DType>& d1, const Vector<DType>& o1,
    const Matrix<DType>& d2, const Vector<DType>& o2)
{
    assert(d1.shape(0) == o1.shape(0));
    assert(d2.shape(0) == o2.shape(0));
    assert(d1.shape(0) == d1.shape(0));
    return qr::solve(hstack(d1, d2 * DType(-1)), o2 - o1);
}

} // namespace splx

} // namespace mxm

#endif // _GEOMETRY_SIMPLEX_H_
