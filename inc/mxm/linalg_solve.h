#if !defined(_LINALG_SOLVE_H)
#define _LINALG_SOLVE_H

#include "mxm/linalg_mat_ref.h"
#include <iostream>
#include <algorithm>


namespace mxm
{

template<typename DType>
Matrix<DType> solveLUTriangle(const Matrix<DType>& mat, const Matrix<DType>& b, bool l_tri)
{
    if(!mat.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    if(mat.shape(0) != b.shape(0))
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    size_t idx_start = l_tri ? 0 : b.shape(0) - 1;
    int step = l_tri ? 1 : - 1;
    Matrix<DType> x(Matrix<DType>::zeros(b.shape()));

    for(size_t w = 0; w < x.shape(1); w++)
    {
        for(size_t i = idx_start, _1 = 0; _1 < b.shape(0); _1 ++, i += step)
        {
            DType off_diag_sum(0);
            for(size_t j = idx_start, _2 = 0; _2 < _1; _2 ++, j += step)
            {
                off_diag_sum += mat(i, j) * x(j,w);
            }
            x(i,w) = (b(i,w) - off_diag_sum) * inv(mat(i,i));
        }
    }

    return x;
}


template<typename DType>
Matrix<DType> solveLowerTriangle(const Matrix<DType>& lower, const Matrix<DType>& b)
{
    return solveLUTriangle(lower, b, 1);
}

template<typename DType>
Matrix<DType> solveUpperTriangle(const Matrix<DType>& upper, const Matrix<DType>& b)
{
    return solveLUTriangle(upper, b, 0);
}

namespace qr
{
// u: column vector, direction
// a: column vectors
template<typename DType>
Matrix<DType> project(const Matrix<DType>& u, const Matrix<DType>& a)
{
    Matrix<DType> u_dir = u.normalized();
    return u_dir.dot(u_dir.T()).dot(a);
}

inline std::vector<std::array<size_t, 2>> upperTrianglizeSequence(size_t cols)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < cols; j++)
        for(size_t i = cols-1; i > j; i--)
            ret.push_back({i,j});
    return ret;
}
inline std::vector<std::array<size_t, 2>> upperHessenbergizeSequence(size_t cols)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < cols; j++)
        for(size_t i = cols-1; i > j+1; i--)
            ret.push_back({i,j});
    return ret;
}
inline std::vector<std::array<size_t, 2>> subdiagonalSequence(const Shape& shape)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < shape[1]-1; j++)
        ret.push_back({j+1,j});
    if(shape[0] > shape[1]) ret.push_back({shape[1], shape[1]-1});
    return ret;
}
inline std::vector<std::array<size_t, 2>> superdiagonalSequence(const Shape& shape)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 1; j < shape[1]; j++)
        ret.push_back({j-1,j});
    return ret;
}
enum TraverseSeq{
    eUpperTrianglize,
    eUpperHessenbergize,
    eSubdiagonal
};

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, Matrix<DType>>::type
givensRotation(const Matrix<DType>& v2)
{
    auto norm_v = v2.norm();
    Matrix<DType> so2({2,2},{
        v2(0,0), v2(1,0),
        -v2(1,0), v2(0,0) });
    so2 *= (decltype(norm_v)(1) / norm_v);
    return so2;
}

template<typename DType, unsigned int N>
Matrix<Hypercomplex<DType, N>>
givensRotation(const Matrix<Hypercomplex<DType, N>>& v2)
{
    auto norm_v = v2.norm();
    Matrix<Hypercomplex<DType, N>> su2({2,2},{
        v2(0,0).conj(), v2(1,0).conj(),
        -v2(1,0), v2(0,0) });
    su2 *= (decltype(norm_v)(1) / norm_v);
    return su2;
}

template<typename DType>
// Matrix<DType> calcMatQFromRotation(const Matrix<DType>& mat_in, bool hessenberg=false)
std::array<Matrix<DType>, 2>
decomposeByRotation(const Matrix<DType>& mat_in, TraverseSeq idx_seq=eUpperTrianglize, bool symmetric=false)
{
    // reference: https://www.math.usm.edu/lambers/mat610/sum10/lecture9.pdf
    if(!mat_in.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    std::array<Matrix<DType>, 2> ret;
    Matrix<DType>& mat(ret[1]);
    mat = mat_in;
    const size_t& n = mat.shape(0);
    Matrix<DType>& rot(ret[0]);
    rot = Matrix<DType>::Identity(n);

    std::vector<std::array<size_t, 2>> seq;
    if(idx_seq == eUpperTrianglize) seq = upperTrianglizeSequence(mat_in.shape(1));
    else if (idx_seq == eUpperHessenbergize) seq = upperHessenbergizeSequence(mat_in.shape(1));
    else if (idx_seq == eSubdiagonal)  seq = subdiagonalSequence(mat_in.shape());

    for(auto& idx: seq)
    {
        auto i = idx[0];
        auto j = idx[1];
        if(norm(mat(i,j)) < eps() * eps()) { continue; }

        Matrix<DType> sub_rot = Matrix<DType>::Identity(n);

        auto so2 = givensRotation(mat(Block({i-1, i+1}, {j, j+1})));

        sub_rot.setBlock(i-1, i-1, so2);

        rot = sub_rot.matmul(rot);
        mat = sub_rot.matmul(mat);
        if(symmetric)
            mat = mat.matmul(sub_rot.T());

        // if(std::is_same<DType, Complex<FloatType>>::value)
        // {
        //     // std::cout << "mat: \n" << mat.str() << std::endl;
        //     std::cout << "recover: \n" << (conj(rot.T())).matmul(mat).str() << std::endl;
        // }
    }
    rot = conj(rot.T());
    return ret;
}

template<typename DType>
Matrix<DType> calcMatQFromReflection(const Matrix<DType>& mat)
{
    //reference: https://en.wikipedia.org/wiki/QR_decomposition
    if(!mat.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    Matrix<DType> mat_u({mat.shape(0), mat.shape(1)});
    for(size_t i = 0; i < mat.shape(0); i++)
    {
        mat_u(Col(i)) = mat(Col(i));
        for(size_t j = 0; j < i; j++)
        {
            mat_u(Col(i)) -= project(mat_u(Col(j)), mat(Col(i)));
        }
    }

    // std::cout << "Bq: \n" << mat_u.str();
    for(size_t i = 0; i < mat.shape(0); i++)
    {
        mat_u(Col(i)) = mat_u(Col(i)).normalized();
    }

    return mat_u;
}

template<typename DType>
Matrix<DType> calcMatQ(const Matrix<DType>& mat)
{
    return decomposeByRotation(mat)[0];
}

template<typename DType>
Matrix<DType> solve(const Matrix<DType>& mat_a, const Matrix<DType>& b)
{
    // Matrix<DType> mat_q(calcMatQ(mat_a));
    auto q_r = decomposeByRotation(mat_a);

    if((q_r[0].matmul(q_r[0])).norm() - q_r[0].shape(0) > eps())
    {
        // singular
        std::cout << "Singular Matrix!" << std::endl;
        return Matrix<DType>::zeros(b.shape());
    }

    // Matrix<DType> mat_r(mat_q.T().matmul(mat_a));
    Matrix<DType> x(solveUpperTriangle(q_r[1], conj(q_r[0].T()).matmul(b)));
    return x;
}

// Reference:
// https://math.stackexchange.com/questions/1262363/convergence-of-qr-algorithm-to-upper-triangular-matrix
template<typename DType>
DType wilkinsonShiftStrategy(const Matrix<DType> mat_a)
{
    size_t n = mat_a.shape(1);
    DType sigma = 0.5 * (mat_a(n-2,n-2) - mat_a(n-1, n-1));
    DType sign = sigma > 0 ? 1. : -1.;
    DType mu = mat_a(n-1, n-1) - (sign * mat_a(n-1, n-2) * mat_a(n-1, n-2)) / (abs(sigma) + sqrt(sigma * sigma) + mat_a(n-1, n-1) * mat_a(n-1, n-1));
    return mu;
}

template<typename DType>
DType errorOrthogonalBlockDiagonal(const Matrix<DType> mat_a)
{
    DType sum = 0.;
    auto seq = superdiagonalSequence(mat_a.shape());
    for(auto & idx : seq)
    {
        DType res = mat_a(idx[0], idx[1]) + mat_a(idx[1], idx[0]);
        sum += res * res;
    }
    return sum;
}

} // namespace qr

template <typename DType>
DType Matrix<DType>::det() const
{
    if(square() && shape(0) == 2)
        return (*this)(0,0) * (*this)(1,1) - (*this)(1,0)*(*this)(0,1);

    Matrix<DType> mat_q(qr::calcMatQ(*this));

    // check full rank
    if((mat_q.matmul(mat_q)).trace() - shape(0) > eps())
        return 0.;

    Matrix<DType> mat_r(mat_q.T().matmul(*this));
    DType det(1);
    for(size_t i = 0; i < shape(0); i++) det *= mat_r(i,i);
    return det;
}

template <typename DType>
Matrix<DType> Matrix<DType>::inv() const
{
    if(!square()) return Matrix<DType>::zeros(shape());
    return qr::solve(*this, Matrix<DType>::Identity(shape(0)));
}

template <typename DType>
std::vector<Complex<DType>> eigvals2x2(const Matrix<DType>& mat)
{
    std::vector<Complex<DType>> ret(2);
    DType tr = mat.trace();
    DType det = mat.det();
    DType delta = tr*tr - 4 * det;
    if(delta >= 0)
    {
        ret.at(0) = Complex<DType>({DType(0.5) * (tr + sqrt(delta)), 0});
        ret.at(1) = Complex<DType>({DType(0.5) * (tr - sqrt(delta)), 0});
    }else
    {
        ret.at(0) = Complex<DType>({DType(0.5) * tr,  DType(0.5) * sqrt(-delta)});
        ret.at(1) = Complex<DType>({DType(0.5) * tr, -DType(0.5) * sqrt(-delta)});
    }
    return ret;
}

// Reference:
// http://www.math.usm.edu/lambers/mat610/sum10/lecture15.pdf
// https://www.cs.purdue.edu/homes/skeel/CS515/4.pdf
template <typename DType>
std::vector<Complex<DType>> eigvals(const Matrix<DType>& mat)
{
    if(!mat.square()) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    size_t n = mat.shape(0);
    DType tol = 1e-13;
    size_t max_it = 30;

    if(2 == n) return eigvals2x2(mat);

    Matrix<DType> hessenberg = qr::decomposeByRotation(mat, qr::eUpperHessenbergize, true)[1];
    Matrix<DType> quasi = hessenberg;

    std::array<Matrix<DType>, 2> q_r;

    for(size_t i = 0; i < max_it; i++)
    {
        // FloatType rho = 1;
        DType rho = qr::wilkinsonShiftStrategy(quasi);
        Matrix<DType> shift = Matrix<DType>::Identity(n) * rho;
        q_r = qr::decomposeByRotation(quasi - shift, qr::eSubdiagonal);
        quasi = q_r[1].matmul(q_r[0]) + shift;
        if(qr::errorOrthogonalBlockDiagonal(q_r[0]) < tol) break;
    }

    // std::cout << quasi.str() << std::endl;

    std::vector<Complex<DType>> ret(n);
    for(size_t i = 0; i < n;)
    {
        if(i < n-1 && abs(q_r[0](i,i) - q_r[0](i+1, i+1)) < tol && abs(q_r[0](i+1,i)) > tol)
        {
            // ret.setBlock(i,0, eigvals2x2(quasi(Block({i,i+2},{i,i+2}))));
            auto eig_pair = eigvals2x2(quasi(Block({i,i+2},{i,i+2})));
            ret.at(i) = eig_pair.at(0);
            ret.at(i+1) = eig_pair.at(1);
            i += 2;
        }else
        {
            ret.at(i) = Complex<DType>({quasi(i,i), 0});
            i++;
        }
    }

    return ret;
}

// References:
// https://en.wikipedia.org/wiki/Inverse_iteration
template<typename DType>
Matrix<DType> inverseIteration(
    const Matrix<DType>& mat_a,
    const DType& eigenvalue,
    const Matrix<DType>& guess)
{
    Matrix<DType> bk = guess;
    size_t n = guess.shape(0);
    typename NormTraits<DType>::type tol = eps() * 20;
    size_t max_it = 20;
    for(size_t i = 0; i < max_it; i++)
    {
        bk = qr::solve(mat_a - eigenvalue * Matrix<DType>::Identity(n), bk);
        bk.normalize();
        auto err = (mat_a.matmul(bk) - eigenvalue * bk).norm();
        // std::cout << "bk: " << bk.T().str() << ", err: " << err << std::endl;
        if(err < tol) break;
    }
    return bk;
}

template<typename DType>
std::array<Matrix<Complex<DType>>, 2>
eig(const Matrix<DType>& mat)
{
    if(!mat.square()) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    size_t n = mat.shape(0);
    std::array<Matrix<Complex<DType>>, 2> ret;
    auto & eig_vals = ret[0];
    auto & eig_vecs = ret[1];
    eig_vals = Matrix<Complex<DType>>({n,1});
    eig_vecs = Matrix<Complex<DType>>({n,n});

    Matrix<Complex<FloatType>> cmat = Matrix<Complex<FloatType>>::zeros({n,n});
    cmat.traverse([&](auto i, auto j){cmat(i,j)(0) = mat(i,j);});

    std::vector<Complex<DType>> prio_que = eigvals(mat);
    std::sort(prio_que.begin(), prio_que.end(), [](auto lhs, auto rhs){return lhs.norm() > rhs.norm();});
    Matrix<Complex<DType>> guess = Matrix<Complex<DType>>::zeros({n,1});
    guess(0,0) = 1;
    for(size_t i = 0;i < prio_que.size(); i++)
    {
        eig_vals(i,0) = prio_que.at(i);
        eig_vecs.setBlock(0,i, inverseIteration(cmat, prio_que.at(i), guess));
    }
    return ret;
}

#if 0
template<typename DType>
std::array<Matrix<Complex<DType>>, 3>
svd(const Matrix<DType>& mat)
{
    std::array<Matrix<Complex<DType>>, 3> ret;
    auto & mat_u = ret[0];
    auto & diag = ret[1];
    auto & mat_vh = ret[2];
    mat_u = Matrix<Complex<DType>>({mat.shape(0), mat.shape(0)});
    diag = Matrix<Complex<DType>>({std::min(mat.shape(0), mat.shape(1))});
    mat_vh = Matrix<Complex<DType>>({n,n});
}
#endif

// References:
// https://www.cs.cornell.edu/courses/cs6210/2010fa/A6/A6.pdf
template <typename DType>
std::array<Matrix<DType>, 2> tridiagonalizeSkewSymmetric(const Matrix<DType>& skew)
{
    assert(("Skew symmetric matrix must be square!", skew.square()));
#if 0 // by reflection
    const size_t& n = skew.shape(0);
    mat_orth = Matrix<DType>::Identity(n);
    for(size_t i = 0; i < n-1; i++)
    {
        Matrix<DType> v = skew(Block({i+1, end()}, {i, i+1}));
        v(0,0) += v.norm();
        v.normalize();
        Matrix<DType> mat_p = Matrix<DType>::Identity(n);
        // std::cout << mat_p(Block({i+1, end()}, {i+1, end()})).shape(0) << "," << v.shape(0) <<  std::endl;
        mat_p(Block({i+1, end()}, {i+1, end()})) -= 2* v.matmul(v.T());
        mat_orth = mat_p.matmul(mat_orth);
        skew = mat_p.matmul(skew).matmul(mat_p);
    }
#else // by rotation

    auto ret = qr::decomposeByRotation(skew, qr::eUpperHessenbergize, /*symmetric*/true);
    return ret;
#endif
}

template <typename DType>
std::array<Matrix<DType>, 2> blockDiagonalizeSkewSymmetric(const Matrix<DType>& skew)
{
    if(!skew.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    std::array<Matrix<DType>, 2> ret;
    Matrix<DType>& mat(ret[1]);
    mat = skew;
    const size_t& n = mat.shape(0);
    Matrix<DType>& rot(ret[0]);
    rot = Matrix<DType>::Identity(n);

    for(size_t j = 0; j < n; j++)
    {
        for(size_t i = n-1; i > (0 == j % 2 ? j+1 : j); i--)
        {
            std::cout << "(i,j): " << i << "," << j << std::endl;
            DType theta = -atan2(mat(i,j), mat(i-1,j));
            // std::cout << "theta: " << theta << std::endl;
            Matrix<DType> sub_rot = Matrix<DType>::Identity(n);

            Matrix<DType> so2({2,2},
                {cos(theta), -sin(theta),
                sin(theta), cos(theta)});

            sub_rot.setBlock(i-1, i-1, so2);
            // std::cout << sub_rot.str() << std::endl;
            std::cout << "b: \n"<< mat.str()  << std::endl;

            rot = sub_rot.matmul(rot);
            mat = sub_rot.matmul(mat).matmul(sub_rot.T());
        }
    }
    rot = rot.T();
    return ret;
}

} // namespace mxm



#endif // _LINALG_SOLVE_H
