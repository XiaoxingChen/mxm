#if !defined(__MODEL_LIDAR_H__)
#define __MODEL_LIDAR_H__

// #include "geometry_ray.h"
#include "rigid_transform.h"
#include "rotation.h"

namespace mxm
{

template<typename DType=float, size_t DIM=3>
class LidarRotating
{
public:
    using ThisType = LidarRotating;

    template<size_t U=DIM, std::enable_if_t<U==3, int> T=0>
    ThisType & setVerticalAngles(const Vector<DType>& angles) { vertical_angles_ = angles; return *this; }

    ThisType & setHorizontalResolution(size_t resolution) { horizontal_resolution_ = resolution; return *this; }
    ThisType & setHorizontalFov(DType start, DType end)
    {
        horizontal_start_ = start;
        horizontal_end_ = end;
        return *this;
    }

    template<size_t U=DIM, std::enable_if_t<U==3, int> T=0>
    Matrix<DType> castRayDirection() const
    {
        Matrix<DType> directions({DIM, vertical_angles_.size() * horizontal_resolution_}, {}, COL);
        DType rot_step = (horizontal_end_ - horizontal_start_) / DType(horizontal_resolution_);
        Matrix<DType> one_shot({DIM, vertical_angles_.size()}, {}, COL);
        for(size_t i = 0; i < vertical_angles_.size(); i++)
        {
            one_shot(0, i) = mxm::cos(vertical_angles_(i));
            one_shot(1, i) = 0;
            one_shot(2, i) = mxm::sin(vertical_angles_(i));
        }

        for(size_t i = 0; i < horizontal_resolution_; i++)
        {
            auto rot = Rotation<DType, DIM>::fromAxisAngle({0,0,1}, horizontal_start_ + DType(i) * rot_step);
            directions.setBlock(0, i * vertical_angles_.size(), rot.apply(one_shot) );
        }
        return pose_.apply(directions);
    }

    template<size_t U=DIM, std::enable_if_t<U==2, int> T=0>
    Matrix<DType> castRayDirection() const
    {
        Matrix<DType> directions({DIM, horizontal_resolution_}, {}, COL);
        DType rot_step = (horizontal_end_ - horizontal_start_) / DType(horizontal_resolution_);
        Matrix<DType> one_shot({DIM,1}, {1., 0.});

        for(size_t i = 0; i < horizontal_resolution_; i++)
        {
            auto rot = Rotation<DType, DIM>::fromAngle(horizontal_start_ + DType(i) * rot_step);
            directions.setBlock(0, i, rot.apply(one_shot) );
        }
        return pose_.apply(directions);
    }

private:
    Vector<DType> vertical_angles_;
    DType horizontal_start_ = (M_PI * -1);
    DType horizontal_end_ = M_PI;
    size_t horizontal_resolution_;
    RigidTransform<DType, DIM> pose_ = RigidTransform<DType, DIM>::identity();
};

} // namespace mxm



#endif // __MODEL_LIDAR_H__
