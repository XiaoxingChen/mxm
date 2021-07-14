#if !defined(__MODEL_CAMERA_H__)
#define __MODEL_CAMERA_H__

#include "geometry_ray.h"
#include "rigid_transform.h"
#include "rotation.h"

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
    virtual Matrix<DType> distort(const Matrix<DType>& pts_3d) const override
    {
        Matrix<DType> ret(pts_3d.shape());
        for(size_t i = 0; i < pts_3d.shape(1); i++)
        {
            assert(norm(pts_3d(2, i) - DType(1)) < std::numeric_limits<DType>::epsilon());

            const auto& x = pts_3d(0, i);
            const auto& y = pts_3d(1, i);
            auto r2 = x*x + y*y;

            auto kr = radial(r2);
            ret(0, i) = x * kr + tangential(0, r2, x, y);
            ret(1, i) = y * kr + tangential(1, r2, x, y);
            ret(2, i) = 1;
        }
        return ret;
    }

    // input points should be in normalized plane
    virtual Matrix<DType> undistort(const Matrix<DType>& pts_3d) const override
    {
        Matrix<DType> ret(pts_3d);
        size_t iter = 5;
        for(size_t i = 0; i < pts_3d.shape(1); i++)
        {
            assert(norm(pts_3d(2, i) - DType(1)) < std::numeric_limits<DType>::epsilon());
            for(size_t j = 0; j < iter; j++)
            {
                const auto& x = pts_3d(0, i);
                const auto& y = pts_3d(1, i);
                auto r2 = x*x + y*y;
                auto kr = radial(r2);
                ret(0, i) = (x - tangential(0, r2, x, y)) / kr;
                ret(1, i) = (y - tangential(1, r2, x, y)) / kr;
                // std::cout << mxm::to_string(ret.T()) << std::endl;
            }

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
    Camera():pose_(RigidTransform<DType>::identity(DIM)), f_(Vec::ones(DIM - 1) * 500.), c_(Vec::ones(DIM - 1) * 300.), resolution_(DIM - 1)
    {
        updateCameraMatrix();
        checkDimension();
    }

    Camera(const RigidTransform<DType>& pose, const Vec& focus=Vec({500, 500}), const Vec& reso=Vec({640, 480}))
        :pose_(pose), f_(focus), c_(reso * 0.5), resolution_(reso)
    {
        updateCameraMatrix();
        checkDimension();
    }
#if 0
    Camera(std::vector<DType>&& focus, std::vector<DType>&& optical_center, std::vector<size_t>&& reso)
        :pose_(RigidTransform<DType>::identity(focus.size() + 1)), f_(std::move(focus)), c_(std::move(optical_center)), resolution_(std::move(reso))
    {
        checkDimension();
        updateCameraMatrix();
    }
#endif

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
    ThisType & setPose(const RigidTransform<DType>& pose) { assert(DIM == pose.dim()); pose_ = pose; return *this; }
    ThisType & setPosition(const Vector<DType>& pos) { assert(DIM == pos.size()); pose_.translation() = pos; return *this; }
    ThisType & setOrientation(const Rotation<DType>& rot) { assert(DIM == rot.dim()); pose_.rotation() = rot; return *this; }
    ThisType & setDistortion(const DistortionPtr<DType>& p) { p_distortion = p; return *this; }

    void operator=(const Camera& rhs)
    {
        pose_ = rhs.pose_;
        f_ = rhs.f_;
        c_ = rhs.c_;
        updateCameraMatrix();
    }

    Ray pixelRay(const std::vector<size_t>& pixel_coordinate) const
    {
        // std::vector<size_t> homo_coord(pixel_coordinate);
        // homo_coord.push_back(1);
        return Ray(pose_.translation(), pixelDirection(Vector<size_t>(pixel_coordinate)));
    }

    Matrix<DType> pixelDirection(const Matrix<size_t>& pixels) const
    {
        Matrix<DType> homo_coord({pixels.shape(0) + 1, pixels.shape(1)});

        pixels.traverse([&](auto i, auto j){homo_coord(i,j) = pixels(i,j);});
        homo_coord(Row(end() - 1)) = Matrix<DType>::ones({1,pixels.shape(1)});
        auto local_dir = cam_mat_inv_.matmul(homo_coord);

        if(p_distortion)
        {
            for(size_t j = 0; j < local_dir.shape(1); j++)
            {
                local_dir(0, j) /= local_dir(2, j);
                local_dir(1, j) /= local_dir(2, j);
                local_dir(2, j) = 1;
            }
            local_dir = p_distortion->distort(local_dir);
        }

        Matrix<DType> directions = pose_.rotation().asMatrix().matmul(local_dir);
        return directions;
    }

    const RigidTransform<DType>& pose() const { return pose_; }
    RigidTransform<DType>& pose() { return pose_; }

#if 0
    std::vector<size_t> resolution() const
    {
        std::vector<size_t> res;
        for(size_t i = 0; i < c_.size(); i++)
            res.push_back( static_cast<size_t>(2 * c_(i) + 0.5) );
        return res;
    }
#endif

    Matrix<DType> project(const Matrix<DType>& points) const
    {
        Matrix<DType> body_frame_points = pose_.inv().apply(points);
        for(size_t i = 0; i < body_frame_points.shape(1); i++)
        {
            body_frame_points(0,i) /= body_frame_points(2,i);
            body_frame_points(1,i) /= body_frame_points(2,i);
            body_frame_points(2,i) = 1;
        }
        if(p_distortion)
        {
            body_frame_points = p_distortion->distort(body_frame_points);
        }

        return cam_mat_.matmul(body_frame_points);
    }

    template<typename PType>
    Matrix<PType> distort(const Matrix<PType>& img_src, bool forward=true)
    {
        if(nullptr == p_distortion) return img_src;
        Matrix<PType> img_out(img_src.shape());
        img_out.traverse([&](auto i, auto j){
            auto dir = cam_mat_inv_.matmul(Vector<DType>{DType(i),DType(j),1});
            auto coord = cam_mat_.matmul(forward ? p_distortion->undistort(dir) : p_distortion->distort(dir));
            img_out(i,j) = interp::bilinearUnitSquare(coord, img_src, "zero")(0,0);
        });
        return img_out;
    }

    // const Matrix<DType>& matrix() const { return cam_mat_; }

    size_t resolution(size_t i) const { return resolution_(i); }

    // Field of View
    // reference:
    // http://kmp.pentaxians.eu/technology/fov/#:~:text=Rectilinear%20Lenses%20on%20Film%20Bodies,the%20diagonal%20of%20the%20film.
    DType fov(size_t axis) const { return 2 * atan2( DType(resolution_(axis)) , 2 * f_(axis)); }
    DType diagFov() const { return 2 * atan2(Vector<DType>(resolution_).norm(), 2 * f_(0)); }

private:

    void updateCameraMatrix()
    {
        cam_mat_ = Matrix<DType>::zeros({pose_.dim() - 1, pose_.dim()});
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
    RigidTransform<DType> pose_ = RigidTransform<DType>::identity(DIM);
    Vector<DType> f_ = Vector<DType>::ones(DIM) * 500.;
    Vector<DType> c_ = Vector<DType>::ones(DIM) * 200.;
    Vector<size_t> resolution_ = Vector<size_t>::ones(DIM) * 400;
    Matrix<DType> cam_mat_;
    Matrix<DType> cam_mat_inv_;
    DistortionPtr<DType> p_distortion = nullptr;
};

} // namespace mxm



#endif // __MODEL_CAMERA_H__
