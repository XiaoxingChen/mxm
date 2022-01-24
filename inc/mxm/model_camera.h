#if !defined(__MODEL_CAMERA_H__)
#define __MODEL_CAMERA_H__

#include "geometry_ray.h"
#include "rigid_transform.h"
#include "rotation.h"
#include "interpolation.h"

#include <memory>

namespace mxm
{
template<typename DType>
class Distortion;
template<typename DType>
using DistortionPtr = std::shared_ptr<Distortion<DType>>;
template<typename DType=float>
class Distortion
{
public:
    virtual Matrix<DType> distort(const Matrix<DType>& pts_3d) const = 0;
    virtual Matrix<DType> undistort(const Matrix<DType>& pts_3d) const = 0;
    static DistortionPtr<DType> radialTangential(const std::vector<DType>& params);
};


// Reference:
// [1] https://docs.opencv.org/3.4/da/d54/group__imgproc__transform.html#ga7dfb72c9cf9780a347fbe3d1c47e5d5a
template<typename DType=float>
class RadialTangentialDistortion: public Distortion<DType>
{
public:
    RadialTangentialDistortion(const std::vector<DType>& v)
    {
        assert(v.size() >= 4);
        k_ = Vector<DType>::zeros(6);
        p_ = Vector<DType>::zeros(2);
        k_(0) = *(v.begin() + 0);
        k_(1) = *(v.begin() + 1);
        p_(0) = *(v.begin() + 2);
        p_(1) = *(v.begin() + 3);
        if(v.size() > 4)
            k_(2) = *(v.begin() + 4);
        if(v.size() > 7)
        {
            k_(3) = *(v.begin() + 5);
            k_(4) = *(v.begin() + 6);
            k_(5) = *(v.begin() + 7);
        }
    }
    DType radial(DType r2) const { return (1 + ((k_(2) * r2 + k_(1))*r2 + k_(0)) * r2) / (1 + ((k_(5) * r2 + k_(4))*r2 + k_(3)) * r2); }
    DType tangential(size_t axis, DType r2, DType x, DType y) const
    {
        if(0 == axis) return 2*p_(0)*x*y + p_(1)*(r2 + 2*x*x);
        if(1 == axis) return 2*p_(1)*x*y + p_(0)*(r2 + 2*y*y);
        return 0;
    }

    // input points should be in normalized plane
    // either {x/z, y/z, 1} or {x/z, y/z}
    virtual Matrix<DType> distort(const Matrix<DType>& homo_pts) const override
    {
        Matrix<DType> ret(homo_pts.shape());
        for(size_t i = 0; i < homo_pts.shape(1); i++)
        {
            assert(2 == homo_pts.shape(0) || norm(homo_pts(2, i) - DType(1)) < std::numeric_limits<DType>::epsilon());

            const auto& x = homo_pts(0, i);
            const auto& y = homo_pts(1, i);
            auto r2 = x*x + y*y;

            auto kr = radial(r2);
            ret(0, i) = x * kr + tangential(0, r2, x, y);
            ret(1, i) = y * kr + tangential(1, r2, x, y);
            if(3 == ret.shape(0)) ret(2, i) = 1;
        }
        return ret;
    }

    // input points should be in normalized plane
    // either {x/z, y/z, 1} or {x/z, y/z}
    virtual Matrix<DType> undistort(const Matrix<DType>& homo_pts) const override
    {
        Matrix<DType> ret(homo_pts);
        const size_t iter_max = 5;
        Matrix<DType> tmp;
        for(size_t i = 0; i < iter_max; i++)
        {
            tmp = distort(ret);
            ret = (homo_pts - tmp) + ret;
        }

        return ret;
    }
private:
    Vector<DType> k_;
    Vector<DType> p_;
};

template<typename DType>
DistortionPtr<DType> Distortion<DType>::radialTangential(const std::vector<DType>& params)
{
    return std::shared_ptr<Distortion<DType>>(new RadialTangentialDistortion<DType>(params));
}
template <typename DType=float, size_t DIM=3>
class Camera
{
public:
    // using DType = float;
    using ThisType = Camera;
    Camera():pose_(RigidTransform<DType, DIM>::identity()), f_(Vec::ones(DIM - 1) * 500.), c_(Vec::ones(DIM - 1) * 300.), resolution_(DIM - 1)
    {
        updateCameraMatrix();
        checkDimension();
    }

    Camera(const RigidTransform<DType, DIM>& pose, const Vec& focus=Vec({500, 500}), const Vec& reso=Vec({640, 480}))
        :pose_(pose), f_(focus), c_(reso * 0.5), resolution_(reso)
    {
        updateCameraMatrix();
        checkDimension();
    }

    ThisType & setFocalLength(const Vector<DType>& f)
    {
        assert(DIM - 1 == f.size());
        f_ = f;
        updateCameraMatrix();
        return *this;
    }
    ThisType & setPrincipalOffset(const Vector<DType>& c) {
        assert(DIM - 1 == c.size());
        c_ = c;
        updateCameraMatrix();
        return *this;
    }
    ThisType & setResolution(const Vector<size_t>& reso) { assert(DIM - 1 == reso.size()); resolution_ = reso; return *this; }
    ThisType & setPose(const RigidTransform<DType, DIM>& pose) { assert(DIM == pose.dim()); pose_ = pose; return *this; }
    ThisType & setPosition(const Vector<DType>& pos) { assert(DIM == pos.size()); pose_.translation() = pos; return *this; }
    ThisType & setOrientation(const Rotation<DType, DIM>& rot) { assert(DIM == rot.dim()); pose_.rotation() = rot; return *this; }
    ThisType & setDistortion(const DistortionPtr<DType>& p) { p_distortion_ = p; return *this; }

    void operator=(const Camera& rhs)
    {
        pose_ = rhs.pose_;
        f_ = rhs.f_;
        c_ = rhs.c_;
        resolution_ = rhs.resolution_;
        p_distortion_ = rhs.p_distortion_;
        updateCameraMatrix();
    }

    Ray pixelRay(const std::vector<size_t>& pixel_coordinate) const
    {
        // std::vector<size_t> homo_coord(pixel_coordinate);
        // homo_coord.push_back(1);
        return Ray(pose_.translation(), pixelDirection(Vector<size_t>(pixel_coordinate)));
    }

    Matrix<DType> pixelDirection(const Matrix<size_t>& pixels, DType z_dir=DType(1.)) const
    {
        auto norm_plane_points =  vstack((pixels - c_) / f_, z_dir * Matrix<DType>::ones({1, pixels.shape(1)}));
        if(p_distortion_)
        {
            norm_plane_points = p_distortion_->undistort(norm_plane_points);
        }
        Matrix<DType> directions = pose_.rotation().apply(norm_plane_points);
        return directions;
    }

    const RigidTransform<DType, DIM>& pose() const { return pose_; }
    RigidTransform<DType, DIM>& pose() { return pose_; }


    Matrix<DType> project(const Matrix<DType>& points) const
    {
        Matrix<DType> body_frame_points = pose_.inv().apply(points);

        Matrix<DType> norm_plane_points = body_frame_points(Block({0, end() - 1}, {})) / body_frame_points(Row(end() - 1));
        if(p_distortion_)
        {
            norm_plane_points = p_distortion_->distort(norm_plane_points);
        }
        return norm_plane_points * f_ + c_;
    }

    template<typename PType>
    Matrix<PType> distort(const Matrix<PType>& img_src, bool forward=true)
    {
        if(nullptr == p_distortion_) return img_src;
        Matrix<PType> img_out(img_src.shape());
        img_out.traverse([&](auto i, auto j){
            auto dir = cam_mat_inv_.matmul(Vector<DType>{DType(i),DType(j),1});
            auto coord = cam_mat_.matmul(forward ? p_distortion_->undistort(dir) : p_distortion_->distort(dir));
            img_out(i,j) = interp::bilinearUnitSquare(coord, img_src, "zero")(0,0);
        });
        return img_out;
    }

    size_t resolution(size_t i) const { return resolution_(i); }

    // Field of View
    // reference:
    // http://kmp.pentaxians.eu/technology/fov/#:~:text=Rectilinear%20Lenses%20on%20Film%20Bodies,the%20diagonal%20of%20the%20film.
    DType fov(size_t axis) const { return 2 * atan2( DType(resolution_(axis)) , 2 * f_(axis)); }
    DType diagFov() const { return 2 * atan2(Vector<DType>(resolution_).norm(), 2 * f_(0)); }

    const Matrix<DType>& invMatrix() const { return cam_mat_inv_; }
    const Matrix<DType>& matrix() const { return cam_mat_; }

private:

    void updateCameraMatrix()
    {
        cam_mat_ = Matrix<DType>::identity(pose_.dim());
        cam_mat_inv_ = Matrix<DType>::identity(pose_.dim());
        for(size_t i = 0; i < pose_.dim() - 1; i++)
        {
            cam_mat_(i,i) = f_(i);
            cam_mat_(i, pose_.dim() - 1) = c_(i);

            cam_mat_inv_(i,i) = 1./ f_(i);
            cam_mat_inv_(i, pose_.dim() - 1) = -c_(i) / f_(i);
        }
    }

    const Camera& checkDimension() const
    {
        assert(DIM == pose_.dim());
        assert(DIM == f_.size() + 1);
        assert(DIM == c_.size() + 1);
        assert(DIM == resolution_.size() + 1);
        assert(DIM == cam_mat_.shape(1));
        assert(DIM == cam_mat_inv_.shape(0));
        assert(DIM == cam_mat_inv_.shape(1));

        return *this;
    }
    RigidTransform<DType, DIM> pose_ = RigidTransform<DType, DIM>::identity();
    Vector<DType> f_ = Vector<DType>::ones(DIM) * 500.;
    Vector<DType> c_ = Vector<DType>::ones(DIM) * 200.;
    Vector<size_t> resolution_ = Vector<size_t>::ones(DIM) * 400;
    Matrix<DType> cam_mat_;
    Matrix<DType> cam_mat_inv_;
    DistortionPtr<DType> p_distortion_ = nullptr;
};

} // namespace mxm



#endif // __MODEL_CAMERA_H__
