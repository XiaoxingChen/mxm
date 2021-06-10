#if !defined(_LINALG_UTILS_H)
#define _LINALG_UTILS_H

#include "linalg_mat.h"
#include <algorithm>

#define NO_MAT_REF 1

namespace mxm
{

template<typename DeriveType>
typename Traits<DeriveType>::DerefType orthogonalComplement(const MatrixBase<DeriveType>& vs_in)
{
    auto & vs = reinterpret_cast<const DeriveType&>(vs_in);
    int diff = static_cast<int>(vs.shape(0)) - vs.shape(1);
    if(diff == -1)
        return orthogonalComplement(vs.T()).T();
    if(diff != 1)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    using DType = typename Traits<DeriveType>::EntryType;
    Matrix<DType> ret({vs.shape(0), 1});
    for(size_t i = 0; i < vs.shape(0); i++)
    {
        Matrix<DType> adjoint_mat(Matrix<DType>::zeros({vs.shape(1), vs.shape(1)}));

        if(i > 0)
            adjoint_mat(Block({0, i},{})) = vs(Block({0, i},{}));
        if(i < vs.shape(0) - 1)
            adjoint_mat(Block({i, },{})) = vs(Block({i+1,},{}));

        ret(i, 0) = mxm::det(adjoint_mat);
        if((i % 2) == 1) ret(i, 0) *= -1;
    }
    return ret;

}

#if 0
inline size_t argMax(const Vec& v)
{
    size_t idx_max(0);
    FloatType val_max(v(0));
    for(size_t i = 0; i < v.size(); i++)
    {
        if(v(i) > val_max) idx_max = i;
    }
    return idx_max;
}
#else
template<typename DType>
typename std::enable_if_t<std::is_arithmetic<DType>::value, std::array<size_t, 2>>
argMax(const Matrix<DType>& mat)
{
    std::array<size_t, 2> ret{0,0};
    DType curr_max = std::numeric_limits<DType>::min();
    mat.traverse([&](auto i, auto j){
        if(mat(i,j) > curr_max)
        {
            curr_max = mat(i,j);
            ret = std::array<size_t, 2>{i,j};
        }
    });
    return ret;
}
#endif

inline std::vector<size_t> argSort(const Vec& v)
{
    std::vector<size_t> indices(v.size());
    sort(indices.begin(), indices.end(),
        [&v](size_t i1, size_t i2){ return v(i1) < v(i2); });
    return indices;
}

inline size_t min(const Vec& v)
{
    FloatType m(INFINITY);
    for(size_t i = 0; i < v.size(); i++) m = v(i) < m ? v(i) : m;
    return m;
}

inline size_t max(const Vec& v)
{
    FloatType m(-INFINITY);
    for(size_t i = 0; i < v.size(); i++) m = v(i) > m ? v(i) : m;
    return m;
}

template <typename DType>
inline Matrix<DType> hstack(const Matrix<DType>& lhs, const Matrix<DType>& rhs)
{
    if(lhs.shape(0) != rhs.shape(0))
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    Matrix<DType> ret({lhs.shape(0), lhs.shape(1) + rhs.shape(1)});
    ret.setBlock(0,0,lhs);
    ret.setBlock(0,lhs.shape(1),rhs);
    return ret;
}

template <typename DType>
inline Matrix<DType> vstack(const Matrix<DType>& lhs, const Matrix<DType>& rhs)
{
    return hstack(lhs.T(), rhs.T()).T();
}

inline size_t factorial(size_t x)
{
    size_t ret = 1;
    while(x > 1) ret *= x--;
    return ret;
}

template<typename DType>
Vector<DType> binaryToVector(size_t dim, uint32_t bin)
{
    Vector<DType> ret(Vector<DType>::zeros(dim));
    size_t axis = 0;
    while(axis < dim)
    {
        if(((bin >> axis) & 1u) > 0) ret(axis) = DType(1);
        axis++;
    }
    return ret;
}

template<typename DType>
Matrix<DType> boundary(const Matrix<DType>& pts)
{
    Matrix<DType> min_max = Matrix<DType>::ones({pts.shape(0), 2}) * INFINITY;
    min_max(Col(1)) *= -1;
    pts.traverse([&](auto i, auto j){
        if(pts(i,j) < min_max(i,0)) min_max(i,0) = pts(i,j);
        if(pts(i,j) > min_max(i,1)) min_max(i,1) = pts(i,j);
    });
    return min_max;
}

template<typename DType>
std::array<DType, 2> elementwiseBounds(const Matrix<DType>& mat)
{
    std::array<DType, 2> min_max{std::numeric_limits<DType>::max(), std::numeric_limits<DType>::min()};
    mat.traverse([&](auto i, auto j){
        min_max[0] = std::min(min_max[0], mat(i,j));
        min_max[1] = std::max(min_max[1], mat(i,j));
    });
    return min_max;
}

template<typename DType>
Matrix<DType> diagonalMatrix(const Matrix<DType>& vec)
{
    Matrix<DType> ret = Matrix<DType>::identity(vec.shape(0));
    for(size_t i = 0; i < vec.shape(0); i++) ret(i,i) = vec(i,0);
    return ret;
}
template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> reduce(const Matrix<DType>& src, const Matrix<CoreType>& core)
{
    Shape ret_shape;
    for(auto i : {0,1}) ret_shape[i] = src.shape(i) / core.shape(i) + size_t(src.shape(i) % core.shape(i));
    Matrix<decltype(DType()*CoreType())> ret(ret_shape);
    ret.traverse([&](auto i, auto j){
        size_t i_src = i * core.shape(0);
        size_t j_src = j * core.shape(1);
        if(i_src + core.shape(0) > src.shape(0) || j_src + core.shape(1) > src.shape(1))
        {
            ret(i,j) = src(i_src, j_src);
            return;
        }
    #if NO_MAT_REF
        ret(i, j) = 0;
        for(size_t u = 0; u < core.shape(0); u++)
        {
            for(size_t v = 0; v < core.shape(1); v++)
            {
                ret(i, j) += core(u,v)* src(i_src + u, j_src + v);
            }
        }
    #else
        ret(i,j) = mxm::sum(src(Block({i_src, i_src + core.shape(0)}, {j_src, j_src + core.shape(1)})) * core);
    #endif
    });
    return ret;
}

template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> convolute(const Matrix<DType>& src, const Matrix<CoreType>& core)
{
    if(core.shape(0) % 2 != 1 || core.shape(1) % 2 != 1)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    Matrix<decltype(DType()*CoreType())> ret(src.shape());

    Shape half_w{core.shape(0)/2, core.shape(1)/2};

    ret.traverse([&](auto i, auto j){
        if(i + half_w[0] >= src.shape(0)
        || j + half_w[1] >= src.shape(1)
        || i < half_w[0]
        || j < half_w[1])
        {
            ret(i,j) = src(i,j);
            return;
        }
#if NO_MAT_REF
        ret(i, j) = 0;
        for(size_t u = 0; u < core.shape(0); u++)
        {
            for(size_t v = 0; v < core.shape(1); v++)
            {
                ret(i, j) += core(u,v)* src(i - half_w[0] + u, j - half_w[1] + v);
            }
        }
#else
        ret(i, j)
            = mxm::sum(core * src(Block(
                {i - half_w[0], i + half_w[0] + 1},
                {j - half_w[1], j + half_w[1] + 1})));
#endif
    });
    return ret;

}
template<typename IntType>
typename std::enable_if_t<std::is_integral<IntType>::value, IntType>
combinations(IntType m, IntType n)
{
    n = std::min(m - n, n);
    IntType denominator = factorial(n);
    IntType numerator = 1;
    for(IntType i = 0; i < n; i++) numerator *= (m - i);
    return numerator / denominator;
}

} // namespace mxm
#endif // _LINALG_UTILS_H
