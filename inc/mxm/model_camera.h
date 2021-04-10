#if !defined(__MODEL_CAMERA_H__)
#define __MODEL_CAMERA_H__

#include "geometry_ray.h"
#include "rigid_transform.h"

namespace mxm
{
class Camera
{
public:
    Camera(size_t dim=3):pose_(RigidTrans::Identity(dim)), f_(Vec::ones(dim - 1) * 500.), c_(Vec::ones(dim - 1) * 300.)
    {
        updateCameraMatrix();
    }
    Camera(const RigidTrans& pose, VecIn focus=Vec({500, 500}), VecIn resolution=Vec({640, 480}))
        :pose_(pose), f_(focus), c_(resolution * 0.5)
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
        if(pose_.dim() > pixel_coordinate.size() + 1)
        {
            std::cout << "position dim: " << pose_.dim()
                << " pixel coordinate dim: " << pixel_coordinate.size() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        Vec dir(pose_.dim());
        for(size_t i = 0; i < pose_.dim() - 1; i++)
        {
            dir(i) = (pixel_coordinate.at(i) - c_(i)) / f_(i);
        }
        dir(pose_.dim() - 1) = 1.;
        dir = pose_.rotation().apply(dir);
        return Ray(pose_.translation(), dir);
    }

    const RigidTrans& pose() const { return pose_; }

    std::vector<size_t> resolution() const
    {
        std::vector<size_t> res;
        for(size_t i = 0; i < c_.size(); i++)
            res.push_back( static_cast<size_t>(2 * c_(i) + 0.5) );
        return res;
    }

    Mat project(const Mat& points) const
    {
        Mat body_frame_points = pose_.inv().apply(points);
        return cam_mat_.matmul(body_frame_points)(Block({0, end() -1}, {}));
    }

private:

    void updateCameraMatrix()
    {
        cam_mat_ = Mat::Identity(pose_.dim());
        for(size_t i = 0; i < pose_.dim() - 1; i++)
        {
            cam_mat_(i,i) = f_(i);
            cam_mat_(i, pose_.dim() - 1) = c_(i);
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
    RigidTrans pose_;
    Vec f_;
    Vec c_;
    Mat cam_mat_;
};

} // namespace mxm



#endif // __MODEL_CAMERA_H__
