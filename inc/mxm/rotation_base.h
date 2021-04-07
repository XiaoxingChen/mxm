#if !defined(_ROTATION_BASE_H_)
#define _ROTATION_BASE_H_

#include "linalg.h"

namespace mxm
{

template<typename DType>
inline Matrix<DType> rodrigues2D(DType angle)
{
    DType c = cos(angle);
    DType s = sin(angle);
    return Matrix<DType>({2,2}, {c, s, -s, c});
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

    return Matrix<DType>::Identity(3) * c + (axis.matmul(axis.T())) * c1 + r_x * s;
}

inline void matrixToAxisAngle3D(const Mat& R, UnitVec& axis, FloatType& angle)
{
    FloatType rx = R(2, 1) - R(1, 2);
    FloatType ry = R(0, 2) - R(2, 0);
    FloatType rz = R(1, 0) - R(0, 1);

    FloatType theta, s, c;
    s = std::sqrt((rx*rx + ry*ry + rz*rz)*0.25);
    c = (R(0, 0) + R(1, 1) + R(2, 2) - 1)*0.5;
    c = c > 1. ? 1. : c < -1. ? -1. : c;
    theta = acos(c);
    Vec r({rx, ry, rz});

    if( s < 1e-5 )
    {
        FloatType t;

        if( c > 0 )
            r = Vec(3);
        else
        {
            t = (R(0, 0) + 1)*0.5;
            r(0) = std::sqrt(std::max(t,(FloatType)0));
            t = (R(1, 1) + 1)*0.5;
            r(1) = std::sqrt(std::max(t,(FloatType)0))*(R(0, 1) < 0 ? -1. : 1.);
            t = (R(2, 2) + 1)*0.5;
            r(2) = std::sqrt(std::max(t,(FloatType)0))*(R(0, 2) < 0 ? -1. : 1.);
            if( fabs(r(0)) < fabs(r(1)) && fabs(r(0)) < fabs(r(2)) && (R(1, 2) > 0) != (r(1)*r(2) > 0) )
                r(2) = -r(2);
            theta /= r.norm();
            r *= theta;
        }
    }
    else
    {
        FloatType vth = 1/(2*s);
        vth *= theta;
        r *= vth;
    }
    angle = r.norm();
    axis = UnitVec(r);
}

inline Mat reflection(UnitVecIn u)
{
    return Mat::Identity(u.size()) - 2 * u.matmul(u.T());
}

inline std::array<UnitVec, 2> planeAngleToBivector(UnitVecIn u, UnitVecIn v, FloatType angle)
{
    UnitVec v_perpend(v - u.dot(v));
    UnitVec v_new(cos(angle/2.) * u + sin(angle/2.) * v_perpend);
    return {u,v_new};
}

inline Mat bivectorToRotationMatrix(UnitVecIn u, UnitVecIn v)
{
    return reflection(v).matmul(reflection(u));
}

inline FloatType normSO2(const Mat& R)
{
    return acos(R(0,0));
}

inline FloatType normSO3(const Mat& R)
{
    UnitVec axis(3);
    FloatType angle;
    matrixToAxisAngle3D(R, axis, angle);
    return angle;
}
} // namespace mxm


#endif // _ROTATION_BASE_H_
