#if !defined(_SPATIAL_AABB_H_)
#define _SPATIAL_AABB_H_

#include "geometry_ray.h"
#include <algorithm>
#include <limits>


namespace mxm
{
class AxisAlignedBoundingBox
{
  public:
    AxisAlignedBoundingBox(Vec min_pt, Vec max_pt):
        min_(min_pt), max_(max_pt) { checkDimension(__FILE__, __LINE__); }
    AxisAlignedBoundingBox(const std::vector<FloatType>& min_pt, const std::vector<FloatType>& max_pt):
        min_(min_pt), max_(max_pt) { checkDimension(__FILE__, __LINE__); }

    AxisAlignedBoundingBox(size_t dim):
        min_(dim), max_(dim) { clear(); }

    using ThisType = AxisAlignedBoundingBox;

    const Vec& min() const {return min_;}
    const Vec& max() const {return max_;}
    Vec center() const {return (min_ + max_) * 0.5;}
    FloatType center(size_t i) const {return (min_(i) + max_(i)) * 0.5;}
    bool empty() const { return min_(0) > max_(0); }
    void clear()
    {
        min_ = Vec::ones(min_.size()) * INFINITY;
        max_ = Vec::ones(max_.size()) * (-INFINITY);
    }
    size_t dim() const { return min_.size(); }
    std::string str() const { return min_.T().str() + max_.T().str(); }

    std::vector<size_t> axesByLength() const { return argSort(max_ - min_); }

    ThisType& extend(const Mat& vertices)
    {
        if(vertices.shape(0) != min_.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        for(size_t i = 0; i < vertices.shape(0); i++)
        {
            for(size_t j = 0; j < vertices.shape(1); j++)
            {
                if(vertices(i, j) > max_(i)) max_(i) = vertices(i, j);
                if(vertices(i, j) < min_(i)) min_(i) = vertices(i, j);
            }
        }
        return *this;
    }

    ThisType& extend(const ThisType& others)
    {
        extend(others.min());
        extend(others.max());
        return *this;
    }

    bool in(const ThisType& rhs) const
    {
        for(size_t i = 0; i < dim(); i++)
        {
            if(min()(i) < rhs.min()(i) || max()(i) > rhs.max()(i)) return false;
        }
        return true;
    }

    void checkDimension(const char* file, uint32_t line) const
    {
        if(min_.size() != max_.size())
            throw std::runtime_error(std::string(file) + ":" + std::to_string(line));
    }

    bool hit(const Ray& ray) const
    {
        if(empty()) return false;
        auto in_out = AxisAlignedBoundingBox::hit(ray, min_, max_);
        // if (in_out[0] < ray.tMin()) return false;
        if (in_out[1] < in_out[0] - eps()) return false;
        if(!ray.valid(in_out[0]) && !ray.valid(in_out[1])) return false;

        // if (min_max[0] < 0 && ray.tMax() < min_max[1]) return false;
        return true;
    }

    static std::array<FloatType, 2> hit(const Ray& ray, const Vec& vertex_min, const Vec& vertex_max)
    {
        if(ray.origin().size() != vertex_min.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(ray.origin().size() != vertex_max.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        FloatType t_in = -INFINITY;
        FloatType t_out = INFINITY;
        for(int i = 0; i < vertex_min.size(); i++)
        {
            if(abs(ray.direction()(i)) < std::numeric_limits<FloatType>::min())
            {
                if(vertex_min(i) < ray.origin()(i) + eps() && ray.origin()(i) < vertex_max(i) + eps())
                    continue;

                t_in = 1;
                t_out = -1;
                break;
            }

            FloatType t0 = (vertex_min(i) - ray.origin()(i)) / ray.direction()(i);
            FloatType t1 = (vertex_max(i) - ray.origin()(i)) / ray.direction()(i);

            t_in = std::max(t_in, std::min(t0, t1));
            t_out = std::min(t_out, std::max(t0, t1));
        }
        return {t_in, t_out};
    }

  private:
    Vec min_;
    Vec max_;
};

using AABB = AxisAlignedBoundingBox;
} // namespace mxm



#endif // _SPATIAL_AABB_H_

