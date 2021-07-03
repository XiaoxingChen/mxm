#if !defined(__MODEL_CAMERA_H__)
#define __MODEL_CAMERA_H__

#include "geometry_ray.h"
#include "rigid_transform.h"

namespace mxm
{
class Camera
{
public:
    using DType = float;
    Camera(size_t dim=3):pose_(RigidTrans::identity(dim)), f_(Vec::ones(dim - 1) * 500.), c_(Vec::ones(dim - 1) * 300.)
    {
        updateCameraMatrix();
    }
    Camera(const RigidTrans& pose, const Vec& focus=Vec({500, 500}), const Vec& reso=Vec({640, 480}))
        :pose_(pose), f_(focus), c_(reso * 0.5)
    {
        checkDimension();
        updateCameraMatrix();
    }

    Camera(std::vector<DType>&& focus, std::vector<DType>&& optical_center, std::vector<size_t>&& reso)
        :pose_(RigidTrans::identity(focus.size() + 1)), f_(std::move(focus)), c_(std::move(optical_center)), resolution_(std::move(reso))
    {
        checkDimension();
        updateCameraMatrix();
    }

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

    Mat pixelDirection(const Matrix<size_t>& pixels) const
    {
        Mat homo_coord({pixels.shape(0) + 1, pixels.shape(1)});

        pixels.traverse([&](auto i, auto j){homo_coord(i,j) = pixels(i,j);});
        homo_coord(Row(end() - 1)) = Mat::ones({1,pixels.shape(1)});

        Mat directions = pose_.rotation().asMatrix().matmul(cam_mat_inv_).matmul(homo_coord);
        return directions;
    }

    const RigidTrans& pose() const { return pose_; }
    RigidTrans& pose() { return pose_; }

#if 0
    std::vector<size_t> resolution() const
    {
        std::vector<size_t> res;
        for(size_t i = 0; i < c_.size(); i++)
            res.push_back( static_cast<size_t>(2 * c_(i) + 0.5) );
        return res;
    }
#endif

    Mat project(const Mat& points) const
    {
        Mat body_frame_points = pose_.inv().apply(points);
        return cam_mat_.matmul(body_frame_points)(Block({0, end() -1}, {}));
    }

    size_t resolution(size_t i) const { return resolution_(i); }
    DType fov(size_t axis) const { return 2 * atan2( DType(resolution_(axis)) , f_(axis)); }

private:

    void updateCameraMatrix()
    {
        cam_mat_ = Mat::identity(pose_.dim());
        cam_mat_inv_ = Mat::identity(pose_.dim());
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

        if(pose_.dim() != f_.size() + 1)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        if(pose_.dim() != c_.size() + 1)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

        return *this;
    }
    RigidTransform<DType> pose_;
    Vector<DType> f_;
    Vector<DType> c_;
    Vector<size_t> resolution_;
    Matrix<DType> cam_mat_;
    Matrix<DType> cam_mat_inv_;
};

} // namespace mxm



#endif // __MODEL_CAMERA_H__
