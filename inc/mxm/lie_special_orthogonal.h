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
bool isValid(const Matrix<DType>& mat, DType tol=std::numeric_limits<DType>::epsilon())
{
    if(! isSquare(mat)) return false;
    return isZero(mat + mat.T(), tol);
}

template<size_t N, typename DType>
std::enable_if_t<3 == N, DType>
findAngle(const Matrix<DType>& skew)
{
    DType x = skew(2,1);
    DType y = skew(0,2);
    DType z = skew(1,0);
    return sqrt(x*x + y*y + z*z);
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

    const DType tol = std::numeric_limits<DType>::epsilon() * 10;
    for(size_t i = 0; i < n;)
    {
        if(i < n-1 && abs(block_diag(i+1,i) - block_diag(i,i+1)) > tol && abs(block_diag(i+1,i) + block_diag(i,i+1)) < tol)
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

template<typename DType>
Vector<DType> vee(const Matrix<DType>& skew)
{
    return Vector<DType>({skew(2,1), skew(0,2), skew(1,0)});
}

template<typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType> wedge(const MatrixBase<DeriveType>& v)
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
    size_t n = static_cast<size_t>(0.5 * (1 + sqrt(1 + 4 * v.shape(0) * 2)) + 0.5);
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
        ret += sin(theta(i)) * skew + (1-cos(theta(i))) * skew.matmul(skew);
    }
    return ret;
}

template <size_t N, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
exp(const Matrix<DType>& skew)
{
    DType theta = findAngle<N>(skew);
    if(abs(theta) < std::numeric_limits<DType>::epsilon())
        return Matrix<DType>::identity(N);

    DType itheta = DType(1) / theta;

    return Matrix<DType>::identity(N) + sin(theta) * itheta * skew + (1 - cos(theta)) * itheta * itheta * skew.matmul(skew);
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

} // namespace so

// Special Orthogonal Lie Group
namespace SO
{

template<typename DType>
bool isValid(const Matrix<DType>& mat, DType tol=std::numeric_limits<DType>::epsilon())
{
    return isIdentity(mat.T().matmul(mat), tol) && norm(mxm::det(mat) - DType(1)) < tol;
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
    return acos(cos_theta);
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
        if(i < n-1 && abs(block_diag(i,i) - block_diag(i+1, i+1)) < tol && abs(block_diag(i+1,i) + block_diag(i,i+1)) < tol)
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
        ret(i,j) = sin(DType(i + 1) * theta(j));
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

template <size_t N, typename DType>
std::enable_if_t<3 == N, Matrix<DType>>
log(const Matrix<DType>& rot_mat)
{
    DType theta = findAngle<N>(rot_mat);
    auto plane = Matrix<DType>({N, 2});

    Matrix<DType> skew = rot_mat - rot_mat.T();
    std::vector<std::pair<size_t, DType>> col_norm(N);
    for(size_t i = 0; i < N; i++) col_norm[i] = {i,skew(Col(i)).norm()};
    std::sort(col_norm.begin(), col_norm.end(), [](auto a, auto b){return a.second > b.second;});
    plane(Col(0)) = skew(Col(col_norm.at(0).first));
    plane(Col(1)) = skew(Col(col_norm.at(1).first));
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
Matrix<DType> slerp(const Matrix<DType>& r0, const Matrix<DType>& r1, DType t)
{
    return pow<3>(r1.matmul(r0.T()), t).matmul(r0);
}

} // namespace SO

} // namespace mxm



#endif // _GROUP_SPECIAL_ORTHOGONAL_H_
