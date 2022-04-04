// References:
// [1] Computing Exponentials of Skew Symmetric Matrices And Logarithms of Orthogonal Matrix
//     https://cs.brynmawr.edu/~dxu/206-2550-2.pdf

#if !defined(_LIE_SPECIAL_ORTHOGONAL_H_)
#define _LIE_SPECIAL_ORTHOGONAL_H_


#include "linalg.h"

#include <limits>

namespace mxm
{

// Special Orthogonal Lie Algebra
namespace so
{
template<typename DType>
bool isValid(const Matrix<DType>& mat, DType* error=nullptr, DType tol=eps<typename Traits<DType>::ArithType>())
{
    if(! isSquare(mat)) return false;
    return isZero(mat + mat.T(), error, tol);
}

template<typename DType>
Matrix<DType> normalized(const Matrix<DType>& input)
{
    return DType(0.5) * (input - input.T());
}

template<size_t N, typename DeriveType>
std::enable_if_t<3 == N, typename Traits<DeriveType>::EntryType>
findAngle(const MatrixBase<DeriveType>& skew)
{
    using DType = typename Traits<DeriveType>::EntryType;
    DType x = skew(2,1);
    DType y = skew(0,2);
    DType z = skew(1,0);
    return mxm::sqrt(x*x + y*y + z*z);
}

template<typename DType>
Vector<DType> findAngles(const Matrix<DType>& skew)
{
    size_t n = skew.shape(0);
    // Vector<DType> ret(n / 2);
    std::vector<DType> ret_data;
    ret_data.reserve(n/2);
    Matrix<DType> orthogonal;
    auto block_diag = realSchurDecomposition(skew, &orthogonal);

    const DType tol = eps<typename Traits<DType>::ArithType>() * 10;
    for(size_t i = 0; i < n;)
    {
        if(i < n-1 && mxm::abs(block_diag(i+1,i) - block_diag(i,i+1)) > tol && mxm::abs(block_diag(i+1,i) + block_diag(i,i+1)) < tol)
        {
            ret_data.push_back(0.5 * (block_diag(i+1,i) - block_diag(i,i+1)));
            i += 2;
        }else
        {
            i++;
        }
    }
    // while(ret_data.size() < n/2) ret_data.push_back(DType(0));
    return Vector<DType>(std::move(ret_data));
}

// convert lower triangle of a skew symmetric matrix to a vector
template <typename DType>
Matrix<DType> unsignedVee(const Matrix<DType>& skew)
{
    std::vector<DType> data;
    for(size_t i = 1; i < skew.shape(0); i++)
    {
        for(size_t j = 0; j < i; j++)
        {
            data.push_back(skew(i,j));
        }
    }
    return Matrix<DType>(fixCol(1), std::move(data));
}

template<typename DeriveType>
Vector<typename Traits<DeriveType>::EntryType>
vee(const MatrixBase<DeriveType>& skew)
{
    return Vector<typename Traits<DeriveType>::EntryType>({skew(2,1), skew(0,2), skew(1,0)});
}

template<typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
wedge(const MatrixBase<DeriveType>& v)
{
    return Matrix<typename Traits<DeriveType>::EntryType>({3,3},
    {0, -v(2,0), v(1,0),
    v(2,0), 0, -v(0,0),
    -v(1,0), v(0,0), 0});
}

// convert lower triangle of a skew symmetric matrix to a vector
template <typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
unsignedWedge(const MatrixBase<DeriveType>& v)
{
    using DType = typename Traits<DeriveType>::EntryType;
    size_t n = static_cast<size_t>(0.5 * (1 + std::sqrt(1 + 4 * v.shape(0) * 2)) + 0.5);
    Matrix<DType> skew({n,n});
    size_t k = 0;
    for(size_t i = 1; i < skew.shape(0); i++)
    {
        for(size_t j = 0; j < i; j++)
        {
            skew(i,j) = v(k, 0);
            skew(j,i) = -v(k, 0);
            k++;
        }
    }
    return skew;
}

// Theorem 2.2
template <typename DType>
Matrix<DType> diannaMatrix(const Vector<DType>& theta)
{
    Matrix<DType> ret({theta.size(), theta.size()});

    ret.traverse([&](auto i, auto j){
        if(0 == i)
        {
            ret(i,j) = theta(j);
        }else
        {
            ret(i,j) = -theta(j) * theta(j) * ret(i - 1,j);
        }
    });

    return ret;
}

template <typename DType>
Matrix<DType> diannaCoeff(const Matrix<DType>& skew)
{
    size_t rot_num = skew.shape(0) / 2;
    size_t dof = (skew.shape(0) - 1) * skew.shape(0) / 2;

    Matrix<DType> ret({rot_num, dof});
    Matrix<DType> skew2 = skew.matmul(skew);
    Matrix<DType> skew_pow = skew;

    for(size_t i = 0; i < rot_num; i++)
    {
        ret(Row(i)) = unsignedVee(skew_pow).T();
        skew_pow = skew2.matmul(skew_pow);
    }
    return ret;
}

// Lemma 2.4
template <typename DType>
Matrix<DType> skewFlattenedGeneralRodrigues(const Matrix<DType>& coeff, const Vector<DType>& theta)
{
    size_t n = coeff.shape(1) / theta.size();
    Matrix<DType> ret = Matrix<DType>::identity(n);
    for(size_t i = 0; i < theta.size(); i++)
    {
        auto skew = unsignedWedge(coeff(Row(i)).T());
        ret += mxm::sin(theta(i)) * skew + (DType(1) - mxm::cos(theta(i))) * skew.matmul(skew);
    }
    return ret;
}

template <size_t N, typename DeriveType>
std::enable_if_t<3 == N, Matrix<typename Traits<DeriveType>::EntryType>>
exp(const MatrixBase<DeriveType>& skew)
{
    using DType = typename Traits<DeriveType>::EntryType;
    using ArithType = typename Traits<DType>::ArithType;
    DType theta = findAngle<N>(skew);
    if(mxm::abs(theta) < eps<typename Traits<DType>::ArithType>())
    // if(mxm::isZero(theta, static_cast<ArithType*>(nullptr), eps<ArithType>()))
        return Matrix<DType>::identity(N) + skew;

    DType itheta = DType(1) / theta;

    return Matrix<DType>::identity(N) + mxm::sin(theta) * itheta * skew + (DType(1.) - mxm::cos(theta)) * itheta * itheta * skew.matmul(skew);
}

template <typename DType>
Matrix<DType>
exp(const Matrix<DType>& skew)
{
    auto angles = findAngles(skew);

    auto dianna_mat = diannaMatrix(angles);
    // std::cout << mxm::to_string(dianna_mat) << std::endl;
    auto dianna_coeff = diannaCoeff(skew);
    // std::cout << mxm::to_string(dianna_coeff) << std::endl;
    auto result =  qr::solve(dianna_mat, dianna_coeff);
    auto rot = skewFlattenedGeneralRodrigues(result, angles);
    return rot;
}

template <typename DType>
Matrix<DType> unsignedWedge(const Matrix<DType>& coeff, const Vector<DType>& theta)
{
    size_t n = coeff.shape(1) / theta.size();
    Matrix<DType> ret = Matrix<DType>::zeros({n,n});
    for(size_t i = 0; i < theta.size(); i++)
    {
        ret += theta(i) * unsignedWedge(coeff(Row(i)).T());
    }
    return ret;
}

// By default left jacob
template <typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
jacob(const MatrixBase<DeriveType>& skew)
{
    using DType = typename Traits<DeriveType>::EntryType;
    static const size_t N(3);
    auto angle = findAngle<N>(skew);
    if(mxm::abs(angle) < eps<typename Traits<DType>::ArithType>()) return Matrix<DType>::identity(N);

    auto i_angle = DType(1.) / angle;
    auto i_angle_2 = i_angle * i_angle;
    Matrix<DType> jac = Matrix<DType>::identity(N)
        + skew * (1 - mxm::cos(angle)) * i_angle_2
        + skew.matmul(skew) * (1 - mxm::sin(angle) * i_angle) * i_angle_2;
    return jac;
}

template <typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
jacobInv(const MatrixBase<DeriveType>& skew)
{
    using DType = typename Traits<DeriveType>::EntryType;
    static const size_t N(3);
    auto angle = findAngle<N>(skew);
    if(mxm::abs(angle) < eps<typename Traits<DType>::ArithType>()) return Matrix<DType>::identity(N);

    auto cot_half = DType(1.) / tan(0.5 * angle);
    auto i_angle = DType(1.) / angle;
    Matrix<DType> j_inv = Matrix<DType>::identity(N) - 0.5 * skew + i_angle * (i_angle - 0.5 * cot_half) * skew.matmul(skew);
    return j_inv;
}

} // namespace so

// Special Orthogonal Lie Group
namespace SO
{

template<typename DType>
bool isValid(const Matrix<DType>& mat, DType* p_error=nullptr, DType tol=std::numeric_limits<DType>::epsilon())
{
    return isIdentity(mat.T().matmul(mat), p_error, tol) && norm(mxm::det(mat) - DType(1)) < tol;
}

template<typename DType>
Matrix<DType> normalized(const Matrix<DType>& input)
{
    auto u_d_vh = svd(input);
    return u_d_vh[0].matmul(u_d_vh[2]);
}

template<typename DType>
Matrix<DType> inv(const Matrix<DType>& mat)
{
    return mat.T();
}

template<size_t N, typename DeriveType>
std::enable_if_t<2 == N, typename Traits<DeriveType>::EntryType>
findAngle(const MatrixBase<DeriveType>& mat)
{
    return atan2(mat(1, 0), mat(0,0));
}

template<size_t N, typename DeriveType>
std::enable_if_t<3 == N, typename Traits<DeriveType>::EntryType>
findAngle(const MatrixBase<DeriveType>& mat)
{
    using DType = typename Traits<DeriveType>::EntryType;
    DType cos_theta = 0.5 * (mat.trace() - N) + 1;
    return mxm::acos(cos_theta);
}

template<typename DType>
Vector<DType> findAngles(const Matrix<DType>& rot)
{
    size_t n = rot.shape(0);
    // Vector<DType> ret(n / 2);
    std::vector<DType> ret_data;
    ret_data.reserve(n/2);
    Matrix<DType> orthogonal;
    auto block_diag = realSchurDecomposition(rot, &orthogonal);

    const DType tol = std::numeric_limits<DType>::epsilon() * 10;

    for(size_t i = 0; i < n;)
    {
        if(i < n-1 && mxm::abs(block_diag(i,i) - block_diag(i+1, i+1)) < tol && mxm::abs(block_diag(i+1,i) + block_diag(i,i+1)) < tol)
        {
            ret_data.push_back(findAngle<2>(block_diag(Block({i,i+2},{i,i+2}))));
            i += 2;
        }else
        {
            i++;
        }
    }
    return Vector<DType>(std::move(ret_data));
}

// Reference:
// chapter 4. formula (5) of paper [1].
template <typename DType>
Matrix<DType> diannaMatrix(const Vector<DType>& theta)
{
    Matrix<DType> ret({theta.size(), theta.size()});
    ret.traverse([&](auto i, auto j){
        ret(i,j) = mxm::sin(DType(i + 1) * theta(j));
    });
    return ret;
}

// Reference:
// chapter 4. formula (5) of paper [1].
template <typename DType>
Matrix<DType> diannaCoeff(const Matrix<DType>& rot)
{
    size_t rot_num = rot.shape(0) / 2;
    size_t dof = (rot.shape(0) - 1) * rot.shape(0) / 2;
    Matrix<DType> ret({rot_num, dof});
    Matrix<DType> rot_pow = rot;
    for(size_t i = 0; i < rot_num; i++)
    {
        auto skew = 0.5 * (rot_pow - rot_pow.T());
        ret(Row(i)) = so::unsignedVee(skew).T();
        rot_pow = rot.matmul(rot_pow);
    }
    return ret;
}

template <size_t N, typename DType>
std::enable_if_t<2 == N, Matrix<DType>>
log(const Matrix<DType>& rot_mat)
{
    auto theta = findAngle<N>(rot_mat);
    return Matrix<DType>({2,2},{0, -theta, theta, 0}, ROW);
}

template<typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
planeOfSingleRotation(const MatrixBase<DeriveType>& rot_mat)
{
    using DType = typename Traits<DeriveType>::EntryType;
    Matrix<DType> skew = rot_mat - rot_mat.T();
    Matrix<DType> col_norm = mxm::sum(skew * skew, 0);
    size_t max_col_idx = argMax(col_norm)[1];

    auto plane = Matrix<DType>({rot_mat.shape(0), 2});

    auto vec_x = skew(Col(max_col_idx)).normalized();
    plane(Col(0)) = vec_x;
    plane(Col(1)) = skew.matmul(vec_x);
    return plane;
}

template <size_t N, typename DeriveType>
std::enable_if_t<3 == N, Matrix<typename Traits<DeriveType>::EntryType>>
log(const MatrixBase<DeriveType>& rot_mat)
{
    using DType = typename Traits<DeriveType>::EntryType;
    DType theta = findAngle<N>(rot_mat);
    if(mxm::abs(theta) < eps<typename Traits<DType>::ArithType>()) return Matrix<DType>::zeros({3,3});

    auto plane = planeOfSingleRotation(rot_mat);
    auto axis = theta * orthogonalComplement(plane).normalized();
    return so::wedge(axis);
}

// Reference:
// chapter 4. of paper [1].
template <typename DType>
Matrix<DType> log(const Matrix<DType>& rot_mat)
{
    auto angles = findAngles(rot_mat);

    auto dianna_mat = diannaMatrix(angles);
    auto dianna_coeff = diannaCoeff(rot_mat);
    auto result = qr::solve(dianna_mat, dianna_coeff);
    auto skew = so::unsignedWedge(result, angles);
    return skew;
}

// denote r = exp(z)
// z = log(r)
// r^y = (e^z)^y = exp(y*z) = exp(y*log(r))
template<size_t N, typename DType>
std::enable_if_t<3 == N,Matrix<DType>>
pow(const Matrix<DType>& r, DType y)
{
    return so::exp<N>(y * SO::log<N>(r));
}

template<typename DType>
Matrix<DType> interp(const Matrix<DType>& r0, const Matrix<DType>& r1, DType t)
{
    return pow<3>(r1.matmul(r0.T()), t).matmul(r0);
}

// left lie derivative:
// d Rp / dR
template<size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
derivPoint(const Matrix<DType>& rot, const Matrix<DType>& p)
{
    return -so::wedge(rot.matmul(p));
}

// left lie derivative:
// d ln(R R_1) / dR
template<size_t N=3, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
derivDistance(const Matrix<DType>& rot, const Matrix<DType>& rot_1)
{
    return so::jacobInv(SO::log<N>(rot.matmul(rot_1)));
}

// template <size_t DIM>
// struct DoF
// {
//     static const size_t val = DIM * (DIM - 1) / 2;
// };
template <size_t DIM>
constexpr size_t dof()
{
    return DIM * (DIM - 1) / 2;
}

} // namespace SO

} // namespace mxm



#endif // _LIE_SPECIAL_ORTHOGONAL_H_
