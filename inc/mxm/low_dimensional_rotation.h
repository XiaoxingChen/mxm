#if !defined(_LOW_DIMENSIONAL_ROTATION_H_)
#define _LOW_DIMENSIONAL_ROTATION_H_

#include "linalg.h"
#include <math.h>
#include <string>
namespace mxm
{
class LowDimensionalRotation
{
public:
    using ThisType = LowDimensionalRotation;

    ~LowDimensionalRotation() {}

    ThisType operator*(const ThisType& rhs) const;

    Mat apply(const Mat& vector) const;

    static ThisType fromMatrix(const Mat& R);
    static ThisType fromAngle(FloatType angle) { return ThisType(angle); }
    static ThisType fromAxisAngle(UnitVecIn axis, FloatType angle) { return ThisType(axis, angle); }
    static ThisType fromPlaneAngle(UnitVecIn u, UnitVecIn v, FloatType angle)
    {
        Mat plane({u.size(),2});
        plane.setBlock(0,0,u);
        plane.setBlock(0,1,v);
        return ThisType(plane, angle);
    }
    Mat asMatrix() const;

    ThisType inv() const { return ThisType(plane_, -angle_); }

    std::string str() const
    {
        return std::string("plane: \n") + plane_.str() + "angle: " + std::to_string(angle_);
    }

    void checkDimension() const
    {
        if(plane_.shape(1) != 2)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    size_t dim() const { return plane_.shape(0); }

    static ThisType Identity(size_t dim)
    {
        Mat plane({dim, 2});
        plane(0,0) = 1;
        plane(1,1) = 1;
        return ThisType(plane, 0.);
    }

private:
    LowDimensionalRotation(FloatType angle)
    :plane_(Mat::Identity(2)), angle_(angle){ }

    LowDimensionalRotation(const Vec& axis, FloatType angle)
    :plane_(orthogonalComplement(axis)), angle_(angle) { checkDimension(); }

    LowDimensionalRotation(const Mat& plane, FloatType angle)
    :plane_(plane), angle_(angle){ checkDimension(); }

    LowDimensionalRotation()
    :plane_(Mat({2, 3},{1,0,0, 0,1,0}).T()), angle_(0){}

    Mat rodrigues() const;
    Mat rodrigues2D() const;
    Mat rodrigues3D() const;

    Mat plane_;
    FloatType angle_;
};


} // namespace mxm

#endif // _LOW_DIMENSIONAL_ROTATION_H_
