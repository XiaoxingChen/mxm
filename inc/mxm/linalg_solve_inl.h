#if !defined(_LINALG_SOLVE_INL_H_)
#define _LINALG_SOLVE_INL_H_

#if defined(MXM_COMPILED_LIB)
#include "linalg_solve.h"
#endif // MXM_COMPILED_LIB

namespace mxm
{
template<typename DType>
Matrix<DType> diagonalMatrix(const Matrix<DType>& vec);
// https://en.wikipedia.org/wiki/Inversion_(discrete_mathematics)
size_t inversionNumber(const std::vector<size_t>& seq)
{
    size_t cnt = 0;
    if(seq.size() < 2) return cnt;
    for(size_t i = 1; i < seq.size(); i++)
    {
        for(size_t j = 0; j < i; j++)
        {
            if(seq.at(j) > seq.at(i)) cnt++;
        }
    }
    return cnt;
}

namespace qr
{

std::vector<std::array<size_t, 2>> upperTrianglizeSequence(size_t cols)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < cols; j++)
        for(size_t i = cols-1; i > j; i--)
            ret.push_back({i,j});
    return ret;
}
std::vector<std::array<size_t, 2>> upperHessenbergizeSequence(size_t cols)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < cols; j++)
        for(size_t i = cols-1; i > j+1; i--)
            ret.push_back({i,j});
    return ret;
}
std::vector<std::array<size_t, 2>> subdiagonalSequence(const Shape& shape)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 0; j < shape[1]-1; j++)
        ret.push_back({j+1,j});
    if(shape[0] > shape[1]) ret.push_back({shape[1], shape[1]-1});
    return ret;
}
std::vector<std::array<size_t, 2>> superdiagonalSequence(const Shape& shape)
{
    std::vector<std::array<size_t, 2>> ret;
    for(size_t j = 1; j < shape[1]; j++)
        ret.push_back({j-1,j});
    return ret;
}

//
// idx_seq: TraverseSeq
// symmetric: output will be QHQ' = mat if enabled.
// output: {Q, R}
// Reference:
// [1] https://www.math.usm.edu/lambers/mat610/sum10/lecture9.pdf
template<typename DeriveType>
std::array<Matrix<typename Traits<DeriveType>::EntryType>, 2>
decomposeByRotation(const MatrixBase<DeriveType>& mat_in, TraverseSeq idx_seq, bool symmetric)
{
    using DType = typename Traits<DeriveType>::EntryType;
    assert(mat_in.square());

    std::array<Matrix<DType>, 2> ret;
    Matrix<DType>& mat(ret[1]);
    mat = mat_in;
    const size_t& n = mat.shape(0);
    Matrix<DType>& rot(ret[0]);
    rot = Matrix<DType>::identity(n);

    std::vector<std::array<size_t, 2>> seq;
    if(idx_seq == eUpperTrianglize) seq = upperTrianglizeSequence(mat_in.shape(1));
    else if (idx_seq == eUpperHessenbergize) seq = upperHessenbergizeSequence(mat_in.shape(1));
    else if (idx_seq == eSubdiagonal)  seq = subdiagonalSequence(mat_in.shape());

    for(auto& idx: seq)
    {
        auto i = idx[0];
        auto j = idx[1];
        if(norm(mat(i,j)) < eps() * eps()) { continue; }

        Matrix<DType> sub_rot = Matrix<DType>::identity(n);

        auto so2 = givensRotation(mat(Block({i-1, i+1}, {j, j+1})));

        sub_rot.setBlock(i-1, i-1, so2);

        rot = sub_rot.matmul(rot);
        mat = sub_rot.matmul(mat);
        if(symmetric)
            mat = mat.matmul(conj(sub_rot.T()));

        // if(std::is_same<DType, Complex<FloatType>>::value)
        // {
        //     // std::cout << "mat: \n" << mxm::to_string(mat) << std::endl;
        //     std::cout << "recover: \n" << mxm::to_string((conj(rot.T())).matmul(mat)) << std::endl;
        // }
    }
    rot = conj(rot.T());
    return ret;
}

template std::array<Matrix<typename Traits<Matrix<double>>::EntryType>, 2> decomposeByRotation(const MatrixBase<Matrix<double>>& mat_in, TraverseSeq idx_seq, bool symmetric);
template std::array<Matrix<typename Traits<Matrix<float>>::EntryType>, 2> decomposeByRotation(const MatrixBase<Matrix<float>>& mat_in, TraverseSeq idx_seq, bool symmetric);
template std::array<Matrix<typename Traits<Matrix<Complex<double>>>::EntryType>, 2> decomposeByRotation(const MatrixBase<Matrix<Complex<double>>>& mat_in, TraverseSeq idx_seq, bool symmetric);
template std::array<Matrix<typename Traits<Matrix<Complex<float>>>::EntryType>, 2> decomposeByRotation(const MatrixBase<Matrix<Complex<float>>>& mat_in, TraverseSeq idx_seq, bool symmetric);

template std::array<Matrix<typename Traits<MatrixRef<double>>::EntryType>, 2> decomposeByRotation(const MatrixBase<MatrixRef<double>>& mat_in, TraverseSeq idx_seq, bool symmetric);
template std::array<Matrix<typename Traits<MatrixRef<float>>::EntryType>, 2> decomposeByRotation(const MatrixBase<MatrixRef<float>>& mat_in, TraverseSeq idx_seq, bool symmetric);
template std::array<Matrix<typename Traits<MatrixRef<Complex<double>>>::EntryType>, 2> decomposeByRotation(const MatrixBase<MatrixRef<Complex<double>>>& mat_in, TraverseSeq idx_seq, bool symmetric);
template std::array<Matrix<typename Traits<MatrixRef<Complex<float>>>::EntryType>, 2> decomposeByRotation(const MatrixBase<MatrixRef<Complex<float>>>& mat_in, TraverseSeq idx_seq, bool symmetric);

//reference: https://en.wikipedia.org/wiki/QR_decomposition
template<typename DType>
Matrix<DType> calcMatQFromReflection(const Matrix<DType>& mat)
{
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

    // std::cout << "Bq: \n" << mxm::to_string(mat_u);
    for(size_t i = 0; i < mat.shape(0); i++)
    {
        mat_u(Col(i)) = mat_u(Col(i)).normalized();
    }

    return mat_u;
}

template Matrix<float> calcMatQFromReflection(const Matrix<float>& mat);
template Matrix<double> calcMatQFromReflection(const Matrix<double>& mat);
template Matrix<Complex<float>> calcMatQFromReflection(const Matrix<Complex<float>>& mat);
template Matrix<Complex<double>> calcMatQFromReflection(const Matrix<Complex<double>>& mat);

template<typename DType>
Matrix<DType> solve(const Matrix<DType>& mat_a, const Matrix<DType>& b)
{
    // Matrix<DType> mat_q(calcMatQ(mat_a));
    if( norm(mat_a) < eps())
    {
        std::cout << "Zero Matrix!" << std::endl;
        return Matrix<DType>::zeros(b.shape());
    }
    auto q_r = decomposeByRotation(mat_a);

    if((q_r[0].matmul(conj(q_r[0].T())) - Matrix<DType>::identity(q_r[0].shape(0)) ).norm()  > 10 * eps())
    {
        // singular
        auto tmp = q_r[0].matmul(conj(q_r[0].T()));
        std::cout << "q.matmul(q.T()): \n" << mxm::to_string(tmp) << std::endl;
        std::cout << "norm: " << (tmp - Matrix<DType>::identity(q_r[0].shape(0)) ).norm() << std::endl;
        std::cout << "Singular Matrix!" << std::endl;
        return Matrix<DType>::zeros(b.shape());
    }

    // Matrix<DType> mat_r(mat_q.T().matmul(mat_a));
    Matrix<DType> x(solveUpperTriangle(q_r[1], conj(q_r[0].T()).matmul(b)));
    return x;
}

template Matrix<float> solve(const Matrix<float>& mat_a, const Matrix<float>& b);
template Matrix<double> solve(const Matrix<double>& mat_a, const Matrix<double>& b);
template Matrix<Complex<float>> solve(const Matrix<Complex<float>>& mat_a, const Matrix<Complex<float>>& b);
template Matrix<Complex<double>> solve(const Matrix<Complex<double>>& mat_a, const Matrix<Complex<double>>& b);


// DType only available for real numbers.
// Reference:
// https://math.stackexchange.com/questions/1262363/convergence-of-qr-algorithm-to-upper-triangular-matrix
template<typename DType>
DType wilkinsonShiftStrategy(const Matrix<DType> mat_a)
{
    size_t n = mat_a.shape(1);
    DType sigma = (mat_a(n-2,n-2) - mat_a(n-1, n-1)) * DType(0.5);
    if(norm(sigma) < eps<typename Traits<DType>::ArithType>()
       && norm(mat_a(n-1, n-1)) < eps<typename Traits<DType>::ArithType>())
            return DType(0);
    DType sign = sigma > 0 ? 1. : -1.;
    DType mu = mat_a(n-1, n-1) - (sign * mat_a(n-1, n-2) * mat_a(n-1, n-2)) / (abs(sigma) + sqrt(sigma * sigma) + mat_a(n-1, n-1) * mat_a(n-1, n-1));
    return mu;
}

template float wilkinsonShiftStrategy(const Matrix<float> mat_a);
template double wilkinsonShiftStrategy(const Matrix<double> mat_a);
// template Complex<float> wilkinsonShiftStrategy(const Matrix<Complex<float>> mat_a);
// template Complex<double> wilkinsonShiftStrategy(const Matrix<Complex<double>> mat_a);

} // namespace qr

template <typename DeriveType>
typename Traits<DeriveType>::EntryType
det(const MatrixBase<DeriveType>& mat)
{
    using DType = typename Traits<DeriveType>::EntryType;
    if(mat.square() && mat.shape(0) == 2)
        return mat(0,0) * mat(1,1) - mat(1,0)*mat(0,1);

    Matrix<DType> mat_q = qr::decomposeByRotation(mat)[0];

    Matrix<DType> mat_r(mat_q.T().matmul(mat));
    DType det(1);
    for(size_t i = 0; i < mat.shape(0); i++) det *= mat_r(i,i);
    return det;
}

template typename Traits<Matrix<float>>::EntryType det(const MatrixBase<Matrix<float>>& mat);
template typename Traits<Matrix<double>>::EntryType det(const MatrixBase<Matrix<double>>& mat);
template typename Traits<Matrix<Complex<float>>>::EntryType det(const MatrixBase<Matrix<Complex<float>>>& mat);
template typename Traits<Matrix<Complex<double>>>::EntryType det(const MatrixBase<Matrix<Complex<double>>>& mat);

template typename Traits<MatrixRef<float>>::EntryType det(const MatrixBase<MatrixRef<float>>& mat);
template typename Traits<MatrixRef<double>>::EntryType det(const MatrixBase<MatrixRef<double>>& mat);
template typename Traits<MatrixRef<Complex<float>>>::EntryType det(const MatrixBase<MatrixRef<Complex<float>>>& mat);
template typename Traits<MatrixRef<Complex<double>>>::EntryType det(const MatrixBase<MatrixRef<Complex<double>>>& mat);

// Inversion
template<typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
inv(const MatrixBase<DeriveType>& mat)
{
    using EntryType = typename Traits<DeriveType>::EntryType;
    using ArithType = typename Traits<DeriveType>::ArithType;
    auto& self = static_cast<const DeriveType&>(mat);
    if(!mat.square()) return Matrix<EntryType>::zeros(self.shape());
    return qr::solve(self, Matrix<EntryType>::identity(self.shape(0)));
}

template Matrix<typename Traits<Matrix<float>>::EntryType> inv(const MatrixBase<Matrix<float>>& mat);
template Matrix<typename Traits<Matrix<double>>::EntryType> inv(const MatrixBase<Matrix<double>>& mat);
template Matrix<typename Traits<Matrix<Complex<float>>>::EntryType> inv(const MatrixBase<Matrix<Complex<float>>>& mat);
template Matrix<typename Traits<Matrix<Complex<double>>>::EntryType> inv(const MatrixBase<Matrix<Complex<double>>>& mat);

template <template <class> class MatrixType, typename DType>
std::array<Complex<typename Traits<DType>::ArithType>,2> eigvals2x2(const MatrixBase<MatrixType<DType>>& mat)
{
    using ArithType = typename Traits<DType>::ArithType;
    std::array<Complex<ArithType>,2> ret;
    DType tr = mat.trace();
    DType det = mxm::det(mat);
    DType delta = tr*tr - 4 * det;
    if(delta >= 0)
    {
        ret.at(0) = Complex<ArithType>({ArithType(0.5) * (tr + sqrt(delta)), 0});
        ret.at(1) = Complex<ArithType>({ArithType(0.5) * (tr - sqrt(delta)), 0});
        // std::cout << "tr: " << tr <<", det: " << det << ", delta: " << delta << std::endl;
        // std::cout << "mat: \n" << mxm::to_string(mat) << std::endl;
    }else
    {
        ret.at(0) = Complex<ArithType>({ArithType(0.5) * tr,  ArithType(0.5) * sqrt(-delta)});
        ret.at(1) = Complex<ArithType>({ArithType(0.5) * tr, -ArithType(0.5) * sqrt(-delta)});
        // std::cout << "<0" << std::endl;
    }
    return ret;
}

template std::array<Complex<float>,2> eigvals2x2(const MatrixBase<Matrix<float>>& mat);
template std::array<Complex<double>,2> eigvals2x2(const MatrixBase<Matrix<double>>& mat);
template std::array<Complex<float>,2> eigvals2x2(const MatrixBase<MatrixRef<float>>& mat);
template std::array<Complex<double>,2> eigvals2x2(const MatrixBase<MatrixRef<double>>& mat);

template<typename DType>
Matrix<DType> shiftedQRIteration(
    const Matrix<DType>& mat,
    Matrix<DType>* p_orthogonal,
    qr::TraverseSeq idx_seq=qr::eUpperTrianglize,
    size_t max_it=40,
    typename Traits<DType>::ArithType tol=eps<typename Traits<DType>::ArithType>())
{
    std::array<Matrix<DType>, 2> q_r;
    Matrix<DType> ret(mat);
    Matrix<DType> eye = Matrix<DType>::identity(mat.shape(0));

    for(size_t i = 0; i < max_it; i++)
    {
        // FloatType rho = 1;
        DType rho = qr::wilkinsonShiftStrategy(ret);
        Matrix<DType> shift = eye * rho;
        q_r = qr::decomposeByRotation(ret - shift, idx_seq);
        ret = q_r[1].matmul(q_r[0]) + shift;
        if(p_orthogonal) *p_orthogonal = p_orthogonal->matmul(q_r[0]);
        // if(qr::errorOrthogonalBlockDiagonal(q_r[0]) < tol) break;
        if(isIdentity(q_r[0], nullptr, tol)) break;
    }
    return ret;
}

template Matrix<float> shiftedQRIteration( const Matrix<float>& mat, Matrix<float>* p_orthogonal, qr::TraverseSeq idx_seq, size_t max_it, float tol);
template Matrix<double> shiftedQRIteration( const Matrix<double>& mat, Matrix<double>* p_orthogonal, qr::TraverseSeq idx_seq, size_t max_it, double tol);
// template Matrix<Complex<float>> shiftedQRIteration( const Matrix<Complex<float>>& mat, Matrix<Complex<float>>* p_orthogonal, qr::TraverseSeq idx_seq, size_t max_it, float tol);
// template Matrix<Complex<double>> shiftedQRIteration( const Matrix<Complex<double>>& mat, Matrix<Complex<double>>* p_orthogonal, qr::TraverseSeq idx_seq, size_t max_it, Complex<double> tol);

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
    typename Traits<DType>::ArithType tol = eps() * 20;
    size_t max_it = 20;
    for(size_t i = 0; i < max_it; i++)
    {
        bk = qr::solve(mat_a - eigenvalue * Matrix<DType>::identity(n), bk);
        bk.normalize();
        auto err = (mat_a.matmul(bk) - eigenvalue * bk).norm();
        // std::cout << "bk: " << mxm::to_string(bk.T()) << ", err: " << err << std::endl;
        if(err < tol) break;
    }
    return bk;
}

template Matrix<float>
inverseIteration<float>(
    const Matrix<float>& mat_a,
    const float& eigenvalue,
    const Matrix<float>& guess);

template Matrix<double>
inverseIteration<double>(
    const Matrix<double>& mat_a,
    const double& eigenvalue,
    const Matrix<double>& guess);

template Matrix<Complex<float>>
inverseIteration<Complex<float>>(
    const Matrix<Complex<float>>& mat_a,
    const Complex<float>& eigenvalue,
    const Matrix<Complex<float>>& guess);

template Matrix<Complex<double>>
inverseIteration<Complex<double>>(
    const Matrix<Complex<double>>& mat_a,
    const Complex<double>& eigenvalue,
    const Matrix<Complex<double>>& guess);

// Reference:
// http://www.math.usm.edu/lambers/mat610/sum10/lecture15.pdf
// https://www.cs.purdue.edu/homes/skeel/CS515/4.pdf
template <typename DType>
std::vector<Complex<DType>> eigvals(const Matrix<DType>& mat)
{
    if(!mat.square()) throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    size_t n = mat.shape(0);
    auto quasi_upper_triangle = realSchurDecomposition(mat);
    // std::cout << mxm::to_string(quasi_upper_triangle, 12) << std::endl;
    std::vector<Complex<DType>> ret(n);
    for(size_t i = 0; i < n;)
    {
        // if(i < n-1 && abs(q_r[0](i,i) - q_r[0](i+1, i+1)) < tol && abs(q_r[0](i+1,i)) > tol)
        if(i < n-1 && abs(quasi_upper_triangle(i+1,i)) > 50 * std::numeric_limits<DType>::epsilon())
        {
            // ret.setBlock(i,0, eigvals2x2(quasi(Block({i,i+2},{i,i+2}))));
            auto eig_pair = eigvals2x2(quasi_upper_triangle(Block({i,i+2},{i,i+2})));
            ret.at(i) = eig_pair.at(0);
            ret.at(i+1) = eig_pair.at(1);
            i += 2;
        }else
        {
            ret.at(i) = Complex<DType>({quasi_upper_triangle(i,i), 0});
            i++;
        }
    }

    return ret;
}

template std::vector<Complex<float>> eigvals(const Matrix<float>& mat);
template std::vector<Complex<double>> eigvals(const Matrix<double>& mat);

// Reference:
// [1] https://dspace.mit.edu/bitstream/handle/1721.1/75282/18-335j-fall-2006/contents/lecture-notes/lec16handout6pp.pdf
// [2] https://www5.in.tum.de/lehre/vorlesungen/konkr_math/WS_11_12/vorl/NumPro_WS1112_Vorlesung_Kapitel_7.pdf#page=18
template <typename DType>
std::array<Matrix<DType>, 2> symmetricEig(const Matrix<DType>& mat)
{
    using ArithType = typename Traits<DType>::ArithType;
    assert(mat.square());// throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    size_t n = mat.shape(0);
    DType tol = 1e-13;
    size_t max_it = 100;

    std::array<Matrix<DType>, 2> val_vec{Matrix<DType>({n,1}), Matrix<DType>({n,n})};

    auto q_r = qr::decomposeByRotation(mat, qr::eUpperHessenbergize, true);
#if 0
    if(std::is_same<DType, Complex<double>>::value)
    {
        std::cout << mxm::to_string(q_r[1], 5) << std::endl;
    }
#endif
    // auto eig_vecs = q_r[0];
    // Matrix<DType> diag = shiftedQRIteration(q_r[1], &eig_vecs, qr::eSubdiagonal, max_it);

    Matrix<ArithType> hessengerg = matrixAtPart(q_r[1], 0); //real matrix

    Matrix<ArithType> eig_vecs_stage2 = Matrix<ArithType>::identity(n);

    Matrix<ArithType> diag = shiftedQRIteration(hessengerg, &eig_vecs_stage2, qr::eSubdiagonal, max_it);
    auto eig_vecs = q_r[0].matmul(eig_vecs_stage2);

    std::vector<size_t> idx_buffer;
    for(size_t i = 0; i < n; i++) idx_buffer.push_back(i);
    std::sort(idx_buffer.begin(), idx_buffer.end(), [&](size_t i, size_t j) { return diag(i,i) > diag(j,j); });

    // std::cout << mxm::to_string(idx_buffer) << std::endl;

    for(size_t i = 0; i < n; i++)
    {
        size_t idx = idx_buffer.at(i);
        val_vec[0](i,0) = diag(idx,idx);
        val_vec[1](Col(i)) = eig_vecs(Col(idx));
    }

    // if(1 == inversionNumber(idx_buffer) % 2) val_vec[1](Col(n-1)) *= -1;
    return val_vec;
}

template std::array<Matrix<float>, 2> symmetricEig(const Matrix<float>& mat);
template std::array<Matrix<double>, 2> symmetricEig(const Matrix<double>& mat);
template std::array<Matrix<Complex<float>>, 2> symmetricEig(const Matrix<Complex<float>>& mat);
template std::array<Matrix<Complex<double>>, 2> symmetricEig(const Matrix<Complex<double>>& mat);

// A = QTQ', where
// Q is a real, orthogonal matrix,
// T is a real, quasi-upper-triangular matrix that has a block upper-triangular structure
// Reference:
// [1] http://www.math.usm.edu/lambers/mat610/sum10/lecture15.pdf
template<typename DType>
Matrix<DType> realSchurDecomposition(const Matrix<DType>& mat, Matrix<DType>* p_orthogonal)
{
    size_t n = mat.shape(0);
    if(n < 3)
    {
        if(p_orthogonal) *p_orthogonal = Matrix<DType>::identity(n);
        return mat;
    }

    auto q_hesson = qr::decomposeByRotation(mat, qr::eUpperHessenbergize, true);
    if(p_orthogonal) *p_orthogonal = q_hesson[0];
    return shiftedQRIteration(q_hesson[1], p_orthogonal, qr::eSubdiagonal);
}

template Matrix<float> realSchurDecomposition(const Matrix<float>& mat, Matrix<float>* p_orthogonal);
template Matrix<double> realSchurDecomposition(const Matrix<double>& mat, Matrix<double>* p_orthogonal);

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

    Matrix<Complex<DType>> cmat(mat);

    std::vector<Complex<DType>> prio_que = eigvals(mat);
    std::sort(prio_que.begin(), prio_que.end(), [](auto lhs, auto rhs){return lhs.norm() > rhs.norm();});
    Matrix<Complex<DType>> guess = Matrix<Complex<DType>>::zeros({n,1});
    guess(0,0) = 1;
    for(size_t i = 0;i < prio_que.size(); i++)
    {
        eig_vals(i,0) = prio_que.at(i);
        eig_vecs.setBlock(0,i, inverseIteration(cmat, prio_que.at(i), guess));
    }
    // std::cout << mxm::to_string(ret) << std::endl;
    return ret;
}

template std::array<Matrix<Complex<float>>, 2> eig(const Matrix<float>& mat);
template std::array<Matrix<Complex<double>>, 2> eig(const Matrix<double>& mat);

template<typename DeriveType>
std::array<Matrix<typename Traits<DeriveType>::EntryType>, 3>
svd(const MatrixBase<DeriveType>& mat)
{
    using DType = typename Traits<DeriveType>::EntryType;
    std::array<Matrix<DType>, 3> u_s_vh;

    size_t n = std::min(mat.shape(0), mat.shape(1));
    u_s_vh[1] = Matrix<DType>({n,1});
    auto inv_sv = Matrix<DType>({n,1});
    bool positive_definite = true;

    auto val_vec = symmetricEig(mat.matmul(mat.T()));
    Matrix<DType> vec_rhs = symmetricEig(mat.T().matmul(mat))[1].T();
    std::vector<size_t> zero_idx;

    u_s_vh[0] = val_vec[1];
    for(size_t i = 0; i < n; i++)
    {
        if(val_vec[0](i,0) < 0 || norm(val_vec[0](i,0)) < std::numeric_limits<typename Traits<DType>::ArithType>::epsilon())
        {
            inv_sv(i,0) = 0;
            positive_definite = false;
            zero_idx.push_back(i);
        }else
        {
            u_s_vh[1](i,0) = sqrt(val_vec[0](i,0));
            inv_sv(i,0) = DType(1.)/u_s_vh[1](i,0);
        }
    }

    auto rough_vh = diagonalMatrix(inv_sv).matmul(u_s_vh[0].T()).matmul(mat);
    for(size_t i = 0; i < rough_vh.shape(0); i++)
    {
        MatrixRef<DType> row_vec = vec_rhs(Row(i));
        if(row_vec.matmul(rough_vh(Row(i)).T())(0,0) < 0)
        {
            row_vec *= -1;
        }
    }
    u_s_vh[2] = vec_rhs;

    return u_s_vh;
}

template std::array<Matrix<typename Traits<Matrix<float>>::EntryType>, 3> svd(const MatrixBase<Matrix<float>>& mat);
template std::array<Matrix<typename Traits<Matrix<double>>::EntryType>, 3> svd(const MatrixBase<Matrix<double>>& mat);
template std::array<Matrix<typename Traits<MatrixRef<float>>::EntryType>, 3> svd(const MatrixBase<MatrixRef<float>>& mat);
template std::array<Matrix<typename Traits<MatrixRef<double>>::EntryType>, 3> svd(const MatrixBase<MatrixRef<double>>& mat);


} // namespace mxm

#endif // _LINALG_SOLVE_INL_H_
