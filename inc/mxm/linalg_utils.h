#if !defined(_LINALG_UTILS_H)
#define _LINALG_UTILS_H

#include "linalg_mat.h"
#include "linalg_solve.h"
#include <algorithm>
#include <thread>

#define NO_MAT_REF 1

namespace mxm
{

template<typename DeriveType>
typename Traits<DeriveType>::DerefType orthogonalComplement(const MatrixBase<DeriveType>& vs_in)
{
    auto & vs = static_cast<const DeriveType&>(vs_in);
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
template<typename DeriveType>
typename std::enable_if_t<
    std::is_arithmetic<
        typename Traits<DeriveType>::EntryType
    >::value,
    std::array<size_t, 2>
>
argMax(const MatrixBase<DeriveType>& mat)
{
    using DType = typename Traits<DeriveType>::EntryType;
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

// template <typename DType>
template<typename DeriveTypeLhs, typename DeriveTypeRhs>
Matrix<typename Traits<DeriveTypeLhs>::EntryType>
hstack(const MatrixBase<DeriveTypeLhs>& lhs, const MatrixBase<DeriveTypeRhs>& rhs)
{
    using DType = typename Traits<DeriveTypeLhs>::EntryType;
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
void reduceBlock(
    const Matrix<DType>& src,
    const Matrix<CoreType>& core,
    Matrix<decltype(DType()*CoreType())>& ret,
    const Shape& start,
    const Shape& shape)
{
    bool major_ax = src.majorAxis();
    bool minor_ax = !major_ax;
    size_t i = 0;
    size_t j = 0;

    for(size_t major_idx = start[major_ax]; major_idx < start[major_ax] + shape[major_ax]; major_idx++)
    {
        for(size_t minor_idx = start[minor_ax]; minor_idx < start[minor_ax] + shape[minor_ax]; minor_idx++)
        {
            if(ROW == major_ax){ i = major_idx; j = minor_idx; }
            else { i = minor_idx; j = major_idx; }

            size_t i_src = i * core.shape(0);
            size_t j_src = j * core.shape(1);
            if(i_src + core.shape(0) > src.shape(0) || j_src + core.shape(1) > src.shape(1))
            {
                ret(i,j) = src(i_src, j_src);
                continue;
            }
            ret(i, j) = 0;
            for(size_t u = 0; u < core.shape(0); u++)
            {
                for(size_t v = 0; v < core.shape(1); v++)
                {
                    ret(i, j) += core(u,v)* src(i_src + u, j_src + v);
                } // v
            } // u
        } // j
    } // i
}

template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> reduce(const Matrix<DType>& src, const Matrix<CoreType>& core)
{
    Shape ret_shape;
    for(auto i : {0,1}) ret_shape[i] = src.shape(i) / core.shape(i) + size_t(src.shape(i) % core.shape(i));
    Matrix<decltype(DType()*CoreType())> ret(ret_shape, {}, src.majorAxis());
    reduceBlock(src, core, ret, {0, 0}, ret_shape);
    return ret;
}

template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> reduceParallel(const Matrix<DType>& src, const Matrix<CoreType>& core)
{
    const size_t THREAD_NUM = 4;
    std::vector<std::thread> threads;
    size_t slice_axis = src.majorAxis() == ROW ? COL : ROW;
    size_t untouched_axis = src.majorAxis();
    Shape ret_shape;
    for(auto i : {0,1}) ret_shape[i] = src.shape(i) / core.shape(i) + size_t(src.shape(i) % core.shape(i));
    Matrix<decltype(DType()*CoreType())> ret(ret_shape, {}, src.majorAxis());

    for(size_t i = 0; i < THREAD_NUM; i++)
    {
        Shape start;
        start[untouched_axis] = 0;
        start[slice_axis] = i * (ret.shape(slice_axis) / THREAD_NUM);

        Shape shape;
        shape[untouched_axis] = ret.shape(untouched_axis);
        shape[slice_axis] = ret.shape(slice_axis) / THREAD_NUM;
        if(i == THREAD_NUM - 1)
            shape[slice_axis] = ret.shape(slice_axis) - start[slice_axis];

        threads.push_back(std::thread([&, st = start, sh = shape](){
            reduceBlock(src, core, ret, st, sh);
        }));
    }
    for(auto & th: threads)
    {
        th.join();
    }
    return ret;
}

template<typename DType, typename CoreType>
void convoluteBlock(
    const Matrix<DType>& src,
    const Matrix<CoreType>& core,
    const Shape& start,
    const Shape& shape,
    Matrix<decltype(DType()*CoreType())>& ret)
{
    bool major_ax = src.majorAxis();
    bool minor_ax = !major_ax;
    Shape half_w{core.shape(0)/2, core.shape(1)/2};
    size_t i = 0;
    size_t j = 0;
    for(size_t major_idx = start[major_ax]; major_idx < start[major_ax] + shape[major_ax]; major_idx++)
    {
        for(size_t minor_idx = start[minor_ax]; minor_idx < start[minor_ax] + shape[minor_ax]; minor_idx++)
        {
            if(ROW == major_ax){ i = major_idx; j = minor_idx; }
            else { i = minor_idx; j = major_idx; }

            if(i + half_w[0] >= src.shape(0)
            || j + half_w[1] >= src.shape(1)
            || i < half_w[0]
            || j < half_w[1])
            {
                ret(i,j) = src(i,j);
                continue;
            }
            ret(i, j) = 0;
            for(size_t u = 0; u < core.shape(0); u++)
            {
                for(size_t v = 0; v < core.shape(1); v++)
                {
                    ret(i, j) += core(u,v)* src(i - half_w[0] + u, j - half_w[1] + v);
                } // v
            } // u
        } // j
    } // i
}

template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> convoluteParallel(const Matrix<DType>& src, const Matrix<CoreType>& core)
{
    const size_t THREAD_NUM = 4;
    std::vector<std::thread> threads;
    size_t slice_axis = src.majorAxis() == ROW ? COL : ROW;
    size_t untouched_axis = src.majorAxis();
    Matrix<decltype(DType()*CoreType())> ret(src.shape(), {}, src.majorAxis());

    for(size_t i = 0; i < THREAD_NUM; i++)
    {
        Shape start;
        start[untouched_axis] = 0;
        start[slice_axis] = i * (src.shape(slice_axis) / THREAD_NUM);

        Shape shape;
        shape[untouched_axis] = src.shape(untouched_axis);
        shape[slice_axis] = src.shape(slice_axis) / THREAD_NUM;
        if(i == THREAD_NUM - 1)
            shape[slice_axis] = src.shape(slice_axis) - start[slice_axis];

        threads.push_back(std::thread([&, st = start, sh = shape](){
            convoluteBlock(src, core, st, sh, ret);
        }));
    }
    for(auto & th: threads)
    {
        th.join();
    }
    return ret;
}

template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> convoluteParallel(const Matrix<DType>& src, const MatrixRef<CoreType>& core)
{
    Matrix<CoreType> deref_core(core);
    return convoluteParallel(src, deref_core);
}

template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> convolute(const Matrix<DType>& src, const Matrix<CoreType>& core)
{
    if(core.shape(0) % 2 != 1 || core.shape(1) % 2 != 1)
        assert(false);
    Matrix<decltype(DType()*CoreType())> ret(src.shape(), {}, src.majorAxis());

    convoluteBlock(src, core, {0,0}, src.shape(), ret);
    return ret;
    // return convoluteParallel(src, core);
}

template<typename DType, typename CoreType>
Matrix<decltype(DType()*CoreType())> convolute(const Matrix<DType>& src, const MatrixRef<CoreType>& core)
{
    Matrix<CoreType> deref_core(core);
    return convolute(src, deref_core);
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

template<typename DeriveType>
typename Traits<DeriveType>::EntryType sum(const MatrixBase<DeriveType>& src)
{
    using EntryType = typename Traits<DeriveType>::EntryType;
    EntryType ret = Traits<EntryType>::zero();
    src.traverse([&](auto i, auto j){
        ret += src(i,j);
    });
    return ret;
}

template<typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
sum(const MatrixBase<DeriveType>& src, size_t axis)
{
    using EntryType = typename Traits<DeriveType>::EntryType;
    if(0 == axis)
    {
        if(1 == src.shape(0)) return src;
        Matrix<EntryType> ret = Matrix<EntryType>::zeros({1, src.shape(1)});
        src.traverse([&](auto i, auto j){ ret(0, j) += src(i,j); });
        return ret;
    }
    // 1 == axis
    return sum(src.T(), 0).T();
}

template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> sin(DType val) { return std::sin(val); }
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> cos(DType val) { return std::cos(val); }
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> asin(DType val) { return std::asin(val); }
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> acos(DType val) { return std::acos(val); }
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> abs(DType val) { return std::abs(val); }
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> sqrt(DType val) { return std::sqrt(val); }
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> pow(DType val, DType exp) { return std::pow(val, exp); }
template<typename DType> std::enable_if_t<std::is_floating_point<DType>::value, DType> exp(DType val) { return std::exp(val); }

template<typename DType>
std::enable_if_t<std::is_floating_point<DType>::value, bool>
isZero(DType val, DType* p_error, DType tol)
{
    DType error = std::abs(val);
    if(p_error) *p_error = error;
    return error < tol;
}



} // namespace mxm
#endif // _LINALG_UTILS_H
