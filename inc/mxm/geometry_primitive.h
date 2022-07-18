#if !defined(_GEOMETRY_PRIMITIVE_H_)
#define _GEOMETRY_PRIMITIVE_H_

#include "linalg.h"
#include "geometry_ray.h"

namespace mxm
{
// return Vec: [t, k1, k2, ..., kn]
// intersect = sum(k) < 1 ? true : false
template<typename DType>
Vector<DType> intersectEquation(const Matrix<DType>& primitive, const Ray<DType>& ray)
{
    assert(primitive.square() && "square marix required!");
    assert(primitive.shape(0) == ray.origin().size() && "input shape mismatch!");

    // (d BA CA) (t k1 k2).T = OA
    // Ax = b

    Vector<DType> b = primitive(Col(0)) - ray.origin();
    Matrix<DType> mat_a(primitive.shape());
    mat_a(Col(0)) = ray.direction();
    for(size_t i = 1; i < primitive.shape(0); i++)
    {
        mat_a(Col(i)) = primitive(Col(0)) - primitive(Col(i));
    }
    return qr::solve(mat_a, b);
}

inline bool validIntersect(const Vec& x)
{
    for(size_t i = 1; i < x.size(); i++)
    {
        if(x(i) < 0 || x(i) > 1) return false;
    }
    // FloatType sum_k = x.norm(1) - fabs(x(0));
    FloatType sum_k = mxm::sum(x) - x(0);
    return sum_k < 1 + eps();
}

//
// triangle: triangle in N-Dimensional space triangle, shape = {N, 3}
// hit_p: hit point in N-Dimensional space triangle, size = N
// triangle_3d: triangle in 2-Dimensional space triangle, shape = {2, 3}
// hit_p_3d: hit point in 2-Dimensional space triangle, size = 2
inline void putTriangleInPlane(
    const Mat& triangle, const Vec& hit_p,
    Mat& triangle_2d, Vec& hit_p_2d)
{
    triangle_2d = Mat({2,3});
    Vec v_ab(triangle(Col(1)) - triangle(Col(0)));
    Vec v_ac(triangle(Col(2)) - triangle(Col(0)));
    Vec dir_ab(v_ab.normalized());
    FloatType l_ab = v_ab.norm();
    // Bx
    triangle_2d(0, 1) = l_ab;
    // Cx, Cy
    triangle_2d(0, 2) = dir_ab.dot(v_ac);
    triangle_2d(1, 2) = (v_ac - triangle_2d(0, 2)* dir_ab).norm();

    hit_p_2d = Vec(2);
    Vec v_ap(hit_p - triangle(Col(0)));
    hit_p_2d(0) = dir_ab.dot(v_ap);
    hit_p_2d(1) = (v_ap - hit_p_2d(0)* dir_ab).norm();
}

template<typename DType>
inline Vector<DType> primitiveNorm(
    const Matrix<DType>& prim, const Ray<DType>& ray)
{
    Matrix<DType> vs = prim(Block({},{1, end()})) - prim(Block({},{0, end() - 1}));
    Vector<DType> norm = orthogonalComplement(vs);
    if(norm.dot(ray.direction()) > 0) norm *= -1;
    return norm.normalized();
}

template<typename DType>
inline Matrix<DType> getPrimitive(
    const Matrix<DType>& vertex_buffer,
    const Matrix<size_t>& vertex_index_buffer,
    size_t primitive_index)
{
    size_t dim = vertex_buffer.shape(0);
    size_t vertex_per_primitive = vertex_index_buffer.shape(0);
    Matrix<DType> ret({dim, vertex_per_primitive});
    for(size_t i = 0; i < vertex_per_primitive; i++)
    {
        size_t vertex_index = vertex_index_buffer(i, primitive_index);
        ret(Col(i)) = vertex_buffer(Col(vertex_index));
    }
    return ret;
}

} // namespace mxm



#endif // _GEOMETRY_PRIMITIVE_H_
