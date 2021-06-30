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

template<typename DType>
inline DType SO2ToAngle(const Matrix<DType>& mat)
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
    std::vector<std::pair<size_t, DType>> col_norm(dim);
    for(size_t i = 0; i < dim; i++) col_norm[i] = {i,skew(Col(i)).norm()};
    std::sort(col_norm.begin(), col_norm.end(), [](auto a, auto b){return a.second > b.second;});
    plane(Col(0)) = skew(Col(col_norm.at(0).first)).normalized();
    plane(Col(1)) = skew(Col(col_norm.at(1).first)).normalized();
    // std::cout << "plane: " << mxm::to_string(plane.T()) << std::endl;
}

template<typename DType>
inline void matrixToAxisAngle3D(const Matrix<DType>& R, Vector<DType>& axis, DType& angle)
{
    Matrix<DType> plane;

    simpleRotationToPlaneAngle(R, plane, angle);
    axis = orthogonalComplement(plane).normalized();
}

inline Mat reflection(const Vec& u_in)
{
    Vec u(u_in.normalized());
    return Mat::identity(u.size()) - 2 * u.matmul(u.T());
}

inline std::array<Vec, 2> planeAngleToBivector(const Vec& u_in, const Vec& v_in, FloatType angle)
{
    Vec u = u_in.normalized();
    Vec v = v_in.normalized();
    Vec v_perpend(v - u.dot(v));
    Vec v_new(cos(angle/2.) * u + sin(angle/2.) * v_perpend);
    return {u,v_new};
}

inline Mat bivectorToRotationMatrix(const Vec& u, const Vec& v)
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
} // namespace so_n




} // namespace mxm


#endif // _ROTATION_BASE_H_
