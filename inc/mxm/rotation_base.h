#if !defined(_ROTATION_BASE_H_)
#define _ROTATION_BASE_H_

#include "linalg.h"
#include <map>

namespace mxm
{

template<typename DType>
inline Matrix<DType> rodrigues2D(DType angle)
{
    DType c = cos(angle);
    DType s = sin(angle);
    return Matrix<DType>({2,2}, {c, -s, s, c}, ROW);
}

template<typename DeriveType>
typename Traits<DeriveType>::EntryType SO2ToAngle(const MatrixBase<DeriveType>& mat)
{
    return atan2(mat(1, 0), mat(0,0));
}

template<typename DType>
inline Matrix<DType> rodrigues3D(const Vector<DType>& axis_in, DType angle)
{
    Vector<DType> axis = axis_in.normalized();
    DType c = cos(angle);
    DType s = sin(angle);
    DType c1 = 1. - c;
    DType itheta = angle ? (1./angle) : 0.;
    const DType & rx = axis(0);
    const DType & ry = axis(1);
    const DType & rz = axis(2);

    Matrix<DType> r_x({3,3},{
        0,  -rz,  ry,
        rz,   0, -rx,
        -ry, rx,   0});

    return Matrix<DType>::identity(3) * c + (axis.matmul(axis.T())) * c1 + r_x * s;
}

//
// Reference:
// 1. https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Log_map_from_SO.283.29_to_.7F.27.22.60UNIQ--postMath-0000000D-QINU.60.22.27.7F.283.29
// 2. https://github.com/opencv/opencv/blob/4.5.2/modules/calib3d/src/calibration.cpp#L388
#if 0
template<typename DType>
inline void matrixToAxisAngle3D(const Matrix<DType>& R, Vector<DType>& axis, DType& angle)
{
    DType rx = R(2, 1) - R(1, 2);
    DType ry = R(0, 2) - R(2, 0);
    DType rz = R(1, 0) - R(0, 1);

    DType theta, s, c;
    s = std::sqrt((rx*rx + ry*ry + rz*rz)*0.25);
    c = (R(0, 0) + R(1, 1) + R(2, 2) - 1)*0.5;
    c = c > 1. ? 1. : c < -1. ? -1. : c;
    theta = acos(c);
    Vector<DType> r({rx, ry, rz});

    if( s < 1e-5 ) // R is very close to identity matrix.
    {
        if( c > 0 ) // trace(R) > 1
        {
            // r = Vec(3);
            axis = Vector<DType>(3);
            angle = 0;
        }
        else
        {
            //recalculate r
            DType t;
            t = (R(0, 0) + 1)*0.5;
            r(0) = std::sqrt(std::max(t,(DType)0));
            t = (R(1, 1) + 1)*0.5;
            r(1) = std::sqrt(std::max(t,(DType)0))*(R(0, 1) < 0 ? -1. : 1.);
            t = (R(2, 2) + 1)*0.5;
            r(2) = std::sqrt(std::max(t,(DType)0))*(R(0, 2) < 0 ? -1. : 1.);
            if( fabs(r(0)) < fabs(r(1)) && fabs(r(0)) < fabs(r(2)) && (R(1, 2) > 0) != (r(1)*r(2) > 0) )
                r(2) = -r(2);
            // theta /= r.norm();
            // r *= theta;
            angle = theta;
            axis = r.normalized();
        }
    }
    else // general cases
    {
        DType vth = 1/(2*s);
        // vth *= theta;
        // r *= vth;
        angle = theta;
        axis = vth * r;
    }
    // angle = r.norm();
    // axis = r.normalized();
}
#endif

template<typename DType>
inline void simpleRotationToPlaneAngle(const Matrix<DType>& R, Matrix<DType>& plane, DType& angle)
{
    size_t dim = R.shape(0);
    plane = Matrix<DType>({dim, 2});

    DType cos_theta = 0.5 * (R.trace() - dim) + 1;
    angle = acos(cos_theta);
    Matrix<DType> skew = R - R.T();
    Matrix<DType> col_norm = mxm::sum(skew * skew, 0);
    size_t max_col_idx = argMax(col_norm)[1];

    auto vec_x = skew(Col(max_col_idx)).normalized();
    plane(Col(0)) = vec_x;
    plane(Col(1)) = skew.matmul(vec_x);
    // std::cout << "plane: " << mxm::to_string(plane.T()) << std::endl;
}

template<typename DType>
void matrixToAxisAngle3D(const Matrix<DType>& R, Vector<DType>& axis, DType& angle)
{
    Matrix<DType> plane;

    simpleRotationToPlaneAngle(R, plane, angle);
    axis = orthogonalComplement(plane).normalized();
}

template<typename DType>
Matrix<DType> reflection(const Vector<DType>& u_in)
{
    Vector<DType> u(u_in.normalized());
    return Matrix<DType>::identity(u.size()) - 2 * u.matmul(u.T());
}

template<typename DType>
inline std::array<Vector<DType>, 2> planeAngleToBivector(const Vector<DType>& u_in, const Vector<DType>& v_in, DType angle)
{
    Vector<DType> u = u_in.normalized();
    Vector<DType> v = v_in.normalized();
    Vector<DType> v_perpend(v - u.dot(v));
    Vector<DType> v_new(cos(angle/2.) * u + sin(angle/2.) * v_perpend);
    return {u,v_new};
}

template<typename DType>
Matrix<DType> bivectorToRotationMatrix(const Vector<DType>& u, const Vector<DType>& v)
{
    return reflection(v).matmul(reflection(u));
}

inline FloatType normSO2(const Mat& R)
{
    return acos(R(0,0));
}

inline FloatType normSO3(const Mat& R)
{
    Vec axis(3);
    FloatType angle;
    matrixToAxisAngle3D(R, axis, angle);
    return angle;
}

template<typename DType>
DType angularPositiveNormalized(DType in)
{
    DType angle = in;
    const DType PERIOD = DType(2 * M_PI);
    while(angle < 0) angle += PERIOD;
    while(angle > PERIOD) angle -= PERIOD;
    return angle;
}

template<typename DType>
DType angularSymmetricNormalized(DType in)
{
    DType angle = in;
    const DType PERIOD = DType(2 * M_PI);
    while(angle < DType(-M_PI)) angle += PERIOD;
    while(angle > DType(M_PI)) angle -= PERIOD;
    return angle;
}

// r * r1 = r2
// r = r2 * r1.T()
template<typename DType>
DType angularDistance(DType theta1, DType theta2)
{
#if 0
    return atan2(sin(theta1 - theta2), cos(theta1 - theta2));
#else
    return angularSymmetricNormalized(theta2 - theta1);

#endif
}
#if 0
namespace so_n
{
template<typename DType>
Matrix<DType> pow(const Matrix<DType>& r, DType y)
{
    assert(r.square() && r.shape(0) <= 3 && "Only SO(2) and SO(3) pow function supported.");
    if(2 == r.shape(0))
    {
        return rodrigues2D(SO2ToAngle(r) * y);
    }

    if(3 == r.shape(0))
    {
        Vector<DType> axis;
        DType angle;
        matrixToAxisAngle3D(r, axis, angle);
        return rodrigues3D(axis, angle * y);
    }

    return Matrix<DType>::identity(r.shape(0));
}

template<typename DType>
Matrix<DType> slerp(const Matrix<DType>& r0, const Matrix<DType>& r1, DType t)
{
    return so_n::pow(r1.matmul(r0.T()), t).matmul(r0);
}

template<typename DType>
Vector<DType> toAngles(const Matrix<DType>& rot)
{
    size_t n = rot.shape(0);
    // Vector<DType> ret(n / 2);
    std::vector<DType> ret_data;
    ret_data.reserve(n/2);
    Matrix<DType> orthogonal;
    auto block_diag = realSchurDecomposition(rot, &orthogonal);

    // std::cout << mxm::to_string(block_diag) << std::endl;

    const DType tol = std::numeric_limits<DType>::epsilon() * 10;

    for(size_t i = 0; i < n;)
    {
        // if(i < n - 1) std::cout << "i: " << i << ", c1: " << abs(block_diag(i,i) - block_diag(i+1, i+1)) << ", c2: " << abs(block_diag(i+1,i) + block_diag(i,i+1)) << std::endl;

        if(i < n-1 && abs(block_diag(i,i) - block_diag(i+1, i+1)) < tol && abs(block_diag(i+1,i) + block_diag(i,i+1)) < tol)
        {
            // std::cout << "get i: " << i << std::endl;
            ret_data.push_back(SO2ToAngle(block_diag(Block({i,i+2},{i,i+2}))));
            i += 2;
        }else
        {
            i++;
        }
    }
    return Vector<DType>(std::move(ret_data));
}

template<typename DType>
Vector<DType> anglesFromSkewSymmetric(const Matrix<DType>& skew)
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

template <typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType> unsignedWedge(const MatrixBase<DeriveType>& v)
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

// References:
// [1] https://cs.brynmawr.edu/~dxu/206-2550-2.pdf
template <typename DType>
Matrix<DType> diannaMatrix(const Vector<DType>& theta)
{
    Matrix<DType> ret({theta.size(), theta.size()});
    ret.traverse([&](auto i, auto j){
        ret(i,j) = sin(DType(i + 1) * theta(j));
    });
    return ret;
}

// Theorem 2.2
template <typename DType>
Matrix<DType> diannaExpMatrix(const Vector<DType>& theta)
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

// References:
// [1] https://cs.brynmawr.edu/~dxu/206-2550-2.pdf
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
        ret(Row(i)) = unsignedVee(skew).T();
        rot_pow = rot.matmul(rot_pow);
    }
    return ret;
}

template <typename DType>
Matrix<DType> diannaExpCoeff(const Matrix<DType>& skew)
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

template <typename DType>
Matrix<DType> skewSymmetricFromCoeff(const Matrix<DType>& coeff, const Vector<DType>& theta)
{
    size_t n = coeff.shape(1) / theta.size();
    Matrix<DType> ret = Matrix<DType>::zeros({n,n});
    for(size_t i = 0; i < theta.size(); i++)
    {
        ret += theta(i) * unsignedWedge(coeff(Row(i)).T());
    }
    return ret;
}

// Lemma 2.4
template <typename DType>
Matrix<DType> rotationFromCoeff(const Matrix<DType>& coeff, const Vector<DType>& theta)
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

template <typename DType>
Matrix<DType> log(const Matrix<DType>& rot_mat)
{
    auto angles = so_n::toAngles(rot_mat);
    // std::cout << mxm::to_string(angles) << std::endl;

    auto dianna_mat = diannaMatrix(angles);
    auto dianna_coeff = diannaCoeff(rot_mat);
    auto result =  qr::solve(dianna_mat, dianna_coeff);
    auto skew = skewSymmetricFromCoeff(result, angles);
    return skew;
}

template <typename DType>
Matrix<DType> exp(const Matrix<DType>& skew)
{
    auto angles = anglesFromSkewSymmetric(skew);

    auto dianna_mat = diannaExpMatrix(angles);
    // std::cout << mxm::to_string(dianna_mat) << std::endl;
    auto dianna_coeff = diannaExpCoeff(skew);
    // std::cout << mxm::to_string(dianna_coeff) << std::endl;
    auto result =  qr::solve(dianna_mat, dianna_coeff);
    auto rot = rotationFromCoeff(result, angles);
    return rot;
}

} // namespace so_n
#endif

template<typename DType>
void toAxisAngle3D(const Quaternion<DType> & q, Vector<DType>& axis, DType& angle)
{
    { //set output identity
        axis = Vector<DType>::zeros(3);
        axis(0) = 1;
        angle = 0;
    }

    DType norm2 = norm(q);
    if(norm2 < std::numeric_limits<DType>::epsilon()) return;

    for(size_t i = 0; i < axis.size(); i++) axis(i) = q(i+1);

    DType im_norm = axis.norm();
    if(im_norm < std::numeric_limits<DType>::epsilon())
    {
        axis = Vector<DType>::zeros(3);
        axis(0) = 1;
        return;
    }

    angle = std::atan2(im_norm, q(0)) * 2;
}

template<typename DType>
Matrix<DType> toSO3(const Quaternion<DType> & q_in)
{
    Vector<DType> axis;
    DType angle;
    toAxisAngle3D(q_in, axis, angle);
    return rodrigues3D(axis, angle);
}

template<typename DType>
Matrix<DType> asRotationMatrix(const Quaternion<DType>& q)
{
    Matrix<DType> lhs(matrixFromLeftQuaternion(q));
    Matrix<DType> rhs(matrixFromRightQuaternion(q.conj()));
    return lhs.matmul(rhs)(Block({1,},{1,}));
}




} // namespace mxm


#endif // _ROTATION_BASE_H_
