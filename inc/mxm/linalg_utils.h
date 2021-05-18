#if !defined(_LINALG_UTILS_H)
#define _LINALG_UTILS_H

#include "mxm/linalg_mat.h"
#include <algorithm>

namespace mxm
{

template<typename DType>
inline Matrix<DType> orthogonalComplement(const Matrix<DType>& vs)
{
    int diff = static_cast<int>(vs.shape(0)) - vs.shape(1);
    if(diff == -1)
        return orthogonalComplement(vs.T()).T();
    if(diff != 1)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    Matrix<DType> ret({vs.shape(0), 1});
    for(size_t i = 0; i < vs.shape(0); i++)
    {
        Matrix<DType> adjoint_mat(Matrix<DType>::zeros({vs.shape(1), vs.shape(1)}));

        if(i > 0)
            adjoint_mat(Block({0, i},{})) = vs(Block({0, i},{}));
        if(i < vs.shape(0) - 1)
            adjoint_mat(Block({i, },{})) = vs(Block({i+1,},{}));

        ret(i, 0) = adjoint_mat.det();
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
    Matrix<DType> ret = Matrix<DType>::Identity(vec.shape(0));
    for(size_t i = 0; i < vec.shape(0); i++) ret(i,i) = vec(i,0);
    return ret;
}

template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> convolute(const Matrix<DType>& src, const Matrix<CoreType>& core)
{
#if 1 // padding
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
        ret(i, j)
            = mxm::sum(core * src(Block(
                {i - half_w[0], i + half_w[0] + 1},
                {j - half_w[1], j + half_w[1] + 1})));
    });
    return ret;

#else
    Matrix<decltype(DType()*CoreType())> ret({src.shape(0) - core.shape(0) + 1, src.shape(1) - core.shape(1) + 1});

    ret.traverse([&](auto i, auto j){
        ret(i,j) = mxm::sum(core * src(Block({i, i + core.shape(0)}, {j, j+ core.shape(1)})));
    });
    return ret;

#endif
}

} // namespace mxm
#endif // _LINALG_UTILS_H
