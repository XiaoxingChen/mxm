#if !defined(_SPATIAL_AABB_H_)
#define _SPATIAL_AABB_H_

#include "geometry_ray.h"
#include <algorithm>
#include <limits>
#include <map>


namespace mxm
{
template<typename DType>
class AxisAlignedBoundingBox
{
  public:
    AxisAlignedBoundingBox<DType>(Vector<DType> min_pt, Vector<DType> max_pt):
        min_(min_pt), max_(max_pt) { checkDimension(__FILE__, __LINE__); }
    // AxisAlignedBoundingBox(const std::vector<DType>& min_pt, const std::vector<DType>& max_pt):
    //     min_(min_pt), max_(max_pt) { checkDimension(__FILE__, __LINE__); }

    AxisAlignedBoundingBox<DType>(size_t dim):
        min_(dim), max_(dim) { clear(); }

    using ThisType = AxisAlignedBoundingBox<DType>;

    const Vec& min() const {return min_;}
    const Vec& max() const {return max_;}
    Vector<DType> center() const {return (min_ + max_) * 0.5;}
    DType center(size_t i) const {return (min_(i) + max_(i)) * 0.5;}
    bool empty() const { return min_(0) > max_(0); }
    void clear()
    {
        min_ = Vec::ones(min_.size()) * INFINITY;
        max_ = Vec::ones(max_.size()) * (-INFINITY);
    }
    size_t dim() const { return min_.size(); }
    std::string str() const { return mxm::to_string(min_.T()) + mxm::to_string(max_.T()); }

    std::vector<size_t> axesByLength() const { return argSort(max_ - min_); }

    ThisType& extend(const Mat& vertices)
    {
        if(vertices.shape(0) != min_.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        Mat bounds = boundary(vertices);

        for(size_t i = 0; i < bounds.shape(0); i++)
        {
            if(bounds(i, 0) < min_(i)) min_(i) = bounds(i, 0);
            if(bounds(i, 1) > max_(i)) max_(i) = bounds(i, 1);
        }
        return *this;
    }

    ThisType& extend(const ThisType& others)
    {
        if(others.empty()) return *this;
        extend(others.min());
        extend(others.max());
        return *this;
    }

    bool in(const ThisType& rhs) const
    {
        if(empty()) return true;
        for(size_t i = 0; i < dim(); i++)
        {
            if(min()(i) < rhs.min()(i) || max()(i) > rhs.max()(i)) return false;
        }
        return true;
    }

    bool contains(const Mat& pts) const
    {
        if(empty()) return false;
        for(size_t i = 0; i < pts.shape(0); i++)
            for(size_t j = 0; j < pts.shape(1); j++)
                if(pts(i,j) < min()(i) || pts(i,j) > max()(i)) return false;
        return true;
    }

    void checkDimension(const char* file, uint32_t line) const
    {
        assert(min_.size() == max_.size());
    }

    bool hit(const Ray<DType>& ray) const
    {
        if(empty()) return false;
        auto in_out = AxisAlignedBoundingBox<DType>::hit(ray, min_, max_);
        // if (in_out[0] < ray.tMin()) return false;
        if (in_out[1] < in_out[0] - eps()) return false;
        if(!ray.valid(in_out[0]) && !ray.valid(in_out[1])) return false;

        // if (min_max[0] < 0 && ray.tMax() < min_max[1]) return false;
        return true;
    }

    static std::array<DType, 2> hit(const Ray<DType>& ray, const Vec& vertex_min, const Vec& vertex_max)
    {
        if(ray.origin().size() != vertex_min.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(ray.origin().size() != vertex_max.size())
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        DType t_in = -INFINITY;
        DType t_out = INFINITY;
        for(int i = 0; i < vertex_min.size(); i++)
        {
            if(abs(ray.direction()(i)) < std::numeric_limits<DType>::min())
            {
                if(vertex_min(i) < ray.origin()(i) + eps() && ray.origin()(i) < vertex_max(i) + eps())
                    continue;

                t_in = 1;
                t_out = -1;
                break;
            }

            DType t0 = (vertex_min(i) - ray.origin()(i)) / ray.direction()(i);
            DType t1 = (vertex_max(i) - ray.origin()(i)) / ray.direction()(i);

            t_in = std::max(t_in, std::min(t0, t1));
            t_out = std::min(t_out, std::max(t0, t1));
        }
        return {t_in, t_out};
    }

  private:
    Vector<DType> min_;
    Vector<DType> max_;
};

template<typename DType>
std::array<DType, 2> distance(const AxisAlignedBoundingBox<DType>& bbox, const Vec& pt)
{
    // AxisAlignedBoundingBox ret(bbox.dim());
    std::array<DType, 2> ret{INFINITY, -INFINITY};

    // DType max_dist = 0;
    for(uint32_t i = 0; i < (1u << bbox.dim()); i++)
    {
        Vector<DType> t = binaryToVector<DType>(bbox.dim(), i);
        Vector<DType> vertex = (-t + 1) * bbox.min() + t * bbox.max();

        ret[1] = std::max(ret[1], (vertex - pt).norm());
    }

    if(bbox.contains(pt))
    {
        ret[0] = 0;
        return ret;
    }

    Vector<DType> min_dist(bbox.dim());
    for(size_t axis = 0; axis < bbox.dim(); axis++)
    {
        if(bbox.min()(axis) < pt(axis) && pt(axis) < bbox.max()(axis))
        {
            min_dist(axis) = 0.;
        }else
        {
            min_dist(axis) = std::min(abs(bbox.min()(axis) - pt(axis)), abs(bbox.max()(axis) - pt(axis)));
        }
    }

    // ret.extend(min_dist);
    ret[0] = min_dist.norm();
    return ret;
}

template<typename DType>
using AABB = AxisAlignedBoundingBox<DType>;


} // namespace mxm



#endif // _SPATIAL_AABB_H_

