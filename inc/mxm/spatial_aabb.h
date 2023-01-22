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

    const Vector<DType>& min() const {return min_;}
    DType min(size_t i) const {return min_(i);}

    const Vector<DType>& max() const {return max_;}
    DType max(size_t i) const {return max_(i);}

    Vector<DType> center() const {return (min_ + max_) * 0.5;}
    DType center(size_t i) const {return (min_(i) + max_(i)) * 0.5;}
    bool empty() const { return min_(0) > max_(0); }
    void clear()
    {
        min_ *= 0;
        min_ += std::numeric_limits<DType>::max();
        max_ *= 0;
        max_ -= std::numeric_limits<DType>::max();
    }
    size_t dim() const { return min_.size(); }
    std::string str() const { return mxm::to_string(min_.T()) + mxm::to_string(max_.T()); }

    std::vector<size_t> axesByLength() const { return argSort(max_ - min_); }

    ThisType& extend(const Mat& vertices)
    {
        assert(vertices.shape(0) == min_.size());

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

    Vector<uint8_t> contains(const Matrix<DType>& pts, DType tol=DType(0)) const
    {
        Vector<uint8_t> ret = Vector<uint8_t>::zeros(pts.shape(1));
        if(empty()) return ret;
        
        for(size_t j = 0; j < pts.shape(1); j++)
        {
            ret(j) = true;
            for(size_t i = 0; i < pts.shape(0); i++)
            {
                ret(j) &= (pts(i,j) > min()(i) - tol && pts(i,j) < max()(i) + tol);
            }
        }
            
        return ret;
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

    static std::array<DType, 2> hit(const Ray<DType>& ray, const Vector<DType>& vertex_min, const Vector<DType>& vertex_max)
    {
        assert(ray.origin().size() == vertex_min.size());
        assert(ray.origin().size() == vertex_max.size());

        DType t_in = -std::numeric_limits<DType>::max();
        DType t_out = std::numeric_limits<DType>::max();
        for(int i = 0; i < vertex_min.size(); i++)
        {
            if(abs(ray.direction()(i)) < std::numeric_limits<DType>::min())
            {
                if(vertex_min(i) < ray.origin()(i) + eps() && ray.origin()(i) < vertex_max(i) + eps())
                    continue;

                t_in = DType(1);
                t_out = DType(0);
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

// parameters:
// bbox: DIM dimension
// pts: (DIM, N), N points with DIM dimension
// ret: (2, N)
template<typename DType>
Matrix<DType> distance(const AxisAlignedBoundingBox<DType>& bbox, const Matrix<DType>& pts)
{
    // AxisAlignedBoundingBox ret(bbox.dim());
    // std::array<DType, 2> ret{INFINITY, -INFINITY};
    Matrix<DType> ret = Matrix<DType>::zeros({2, pts.shape(1)});
    ret(Row(0)) += std::numeric_limits<DType>::max();
    ret(Row(1)) -= std::numeric_limits<DType>::max();
    auto is_contained = bbox.contains(pts);

    for(size_t pt_idx = 0; pt_idx < pts.shape(1); pt_idx++)
    {
        for(uint32_t i = 0; i < (1u << bbox.dim()); i++)
        {
            Vector<DType> t = binaryToVector<DType>(bbox.dim(), i);
            Vector<DType> vertex = (-t + 1) * bbox.min() + t * bbox.max();

            ret(1, pt_idx) = std::max(ret(1, pt_idx), (vertex - pts(Col(pt_idx))).norm());
        }

        if(is_contained(pt_idx) > 0)
        {
            ret(0, pt_idx) = 0;
            continue;
        }

        Vector<DType> min_dist(bbox.dim());
        for(size_t axis = 0; axis < bbox.dim(); axis++)
        {
            if(bbox.min(axis) < pts(axis, pt_idx) && pts(axis, pt_idx) < bbox.max(axis))
            {
                min_dist(axis) = 0.;
            }else
            {
                min_dist(axis) = std::min(abs(bbox.min(axis) - pts(axis, pt_idx)), abs(bbox.max(axis) - pts(axis, pt_idx)));
            }
        }
        ret(0, pt_idx) = min_dist.norm();
    }
    
    return ret;
}

template<typename DType>
using AABB = AxisAlignedBoundingBox<DType>;


} // namespace mxm



#endif // _SPATIAL_AABB_H_

