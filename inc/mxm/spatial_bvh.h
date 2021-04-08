#if !defined(__SPATIAL_BVH_H__)
#define __SPATIAL_BVH_H__

#include "spatial_aabb.h"
#include "geometry_primitive.h"
#include <stack>
#include <memory>


namespace mxm
{

namespace bvh
{

enum HitType
{
    eClosestHit = 0,
    eMultiHit = 1,
    eAnyHit = 2,
    eHitTypeNum
};

struct Node
{
    Node(size_t dim): aabb(dim), is_leaf(false){}
    AABB aabb;
    std::vector<size_t> children_index_buffer;
    std::vector<size_t> primitive_index_buffer;
    bool is_leaf;
};

class PrimitiveMeshTree
{
public:
    struct HitRecord
    {
        HitRecord(size_t dim): ray(dim){}
        Ray ray;
        FloatType t;
        size_t prim_idx;
    };

    void build(size_t primitive_per_leaf=4, bool verbose=true);
    size_t dim() const { return vertex_buffer_->shape(0); }

    PrimitiveMeshTree(
        std::shared_ptr<Mat>& vertex_buffer,
        std::shared_ptr<Matrix<size_t>>& vertex_index_buffer)
        :vertex_buffer_(vertex_buffer), vertex_index_buffer_(vertex_index_buffer)
        {}

    Mat primitive(size_t idx) const { return getPrimitive(*vertex_buffer_, *vertex_index_buffer_, idx); }

    size_t multiHit(const Ray& ray) const;
    std::vector<HitRecord> hit(const Ray& ray, HitType hit_type) const;

    // const AABB& aabb() const { return node_buffer_.at(0).aabb; }
    const Mat & vertexBuffer() const { return *vertex_buffer_; }
    const size_t primitiveSize() const { return vertex_index_buffer_->shape(1); }

private:

    std::shared_ptr<Mat> vertex_buffer_;
    std::shared_ptr<Matrix<size_t>> vertex_index_buffer_;
    std::vector<bvh::Node> node_buffer_;

};


} // namespace bvh


} // namespace mxm

#ifdef MXM_HEADER_ONLY
#include "spatial_bvh_inl.h"
#endif



#endif // __SPATIAL_BVH_H__
