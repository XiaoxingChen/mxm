#ifndef __GEOMETRY_RAY_H__
#define __GEOMETRY_RAY_H__

#include "linalg.h"
// #include "common.h"
#include "rigid_transform.h"
#include "transform_affine.h"

namespace mxm{

template<typename DType=float>
class Ray
{
public:
    using ThisType = Ray<DType>;
    Ray():origin_(), direction_(){}
    Ray(size_t dimension):origin_(dimension), direction_(dimension){ direction_(0) = 1; }
    Ray(const Vector<DType>& origin, const Vector<DType>& direction, DType t_min=DType(1e-4), DType t_max=DType(1e4)):
    origin_(origin), direction_(direction.normalized()) { checkDimension(); }
    // Ray(const std::vector<DType>& origin, const std::vector<DType>& direction, DType t_min=mxm::tMin(), DType t_max=mxm::tMax()):
    // origin_(origin), direction_(direction), t_min_(t_min), t_max_(t_max) { checkDimension(); }

    const ThisType& checkDimension(size_t dim=0) const
    {
        if(origin_.size() != direction_.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(dim != 0 && dim != origin_.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        return *this;
    }

    const Vector<DType>& origin() const       { return origin_; }
    Vector<DType>& origin()       { return origin_; }

    const Vector<DType>& direction() const    { return direction_; }
    void setDirection(const Vector<DType>& direction)    { direction_ = direction.normalized(); checkDimension(); }

    Vector<DType> operator() (DType t) const { return origin() + direction() * t; }
    bool valid(DType t) const { return t < t_max_ && t > t_min_; }
    DType tMax() const { return t_max_; }
    DType& tMax() { return t_max_; }
    DType tMin() const { return t_min_; }

    private:
    Vector<DType> origin_;
    Vector<DType> direction_;
    DType t_min_ = DType(1e-4);
    DType t_max_ = DType(1e4);
};

template<typename DType, size_t DIM=3>
Ray<DType> apply(const RigidTransform<DType, DIM>& tf, const Ray<DType>& r)
{
    return Ray<DType>(tf.apply(r.origin()), tf.rotation().apply(r.direction()));
}

template<typename DType, size_t DIM=3>
Ray<DType> apply(const AffineTransform<DType, DIM>& tf, const Ray<DType>& r)
{
    Vector<DType> new_direction = tf.linear().matmul(r.direction());
    DType scale = mxm::norm(new_direction);
    return Ray<DType>(tf.apply(r.origin()), new_direction * (DType(1.) / scale), r.tMin() * scale, r.tMax() * scale);
}

// using RayPtr = std::shared_ptr<Ray>;
}//ray tracing

#endif