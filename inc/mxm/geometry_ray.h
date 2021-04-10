#ifndef __GEOMETRY_RAY_H__
#define __GEOMETRY_RAY_H__

#include "linalg.h"
#include "common.h"
#include "rigid_transform.h"

namespace mxm{

class Ray
{
public:
    Ray():origin_(), direction_(), t_min_(mxm::tMin()), t_max_(mxm::tMax()){}
    Ray(size_t dimension):origin_(dimension), direction_(dimension), t_min_(mxm::tMin()), t_max_(mxm::tMax()){}
    Ray(const Vec& origin, const Vec& direction, FloatType t_min=mxm::tMin(), FloatType t_max=mxm::tMax()):
    origin_(origin), direction_(direction), t_min_(t_min), t_max_(t_max) { checkDimension(); }
    Ray(const std::vector<FloatType>& origin, const std::vector<FloatType>& direction, FloatType t_min=mxm::tMin(), FloatType t_max=mxm::tMax()):
    origin_(origin), direction_(direction), t_min_(t_min), t_max_(t_max) { checkDimension(); }

    const Ray& checkDimension(size_t dim=0) const
    {
        if(origin_.size() != direction_.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(dim != 0 && dim != origin_.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        return *this;
    }

    const Vec& origin() const       { return origin_; }
    Vec& origin()       { return origin_; }

    const Vec& direction() const    { return direction_; }
    void setDirection(const Vec& direction)    { direction_ = direction.normalized(); checkDimension(); }

    Vec operator() (FloatType t) const { return origin() + direction() * t; }
    bool valid(FloatType t) const { return t < t_max_ && t > t_min_; }
    FloatType tMax() const { return t_max_; }
    FloatType& tMax() { return t_max_; }
    FloatType tMin() const { return t_min_; }

    private:
    Vec origin_;
    Vec direction_;
    FloatType t_min_;
    FloatType t_max_;
};

inline Ray apply(const RigidTrans& tf, const Ray& r)
{
    return Ray(tf.apply(r.origin()), tf.rotation().apply(r.direction()));
}

// using RayPtr = std::shared_ptr<Ray>;
}//ray tracing

#endif