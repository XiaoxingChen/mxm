#if !defined(__CV_CALIBRATION_H__)
#define __CV_CALIBRATION_H__

#include "cv_basic.h"
#include "linalg_utils.h"
#include "optimize.h"
#include "optimize_sfm_gn.h"

namespace mxm
{

template<typename DType>
class PinholeCameraIntrinsicEstimator
{

public:
    static const size_t POSE_DOF = 6;
    static const size_t DIM = 3;
    static const size_t INTRINSIC_PARAM_NUM = 8; //fx fy cx cy k1 k2 p1 p2

    PinholeCameraIntrinsicEstimator(
        const Matrix<DType>& pts3d,
        const Vector<Matrix<DType>>& pts2d)
    :pts3d_(pts3d), pts2d_(pts2d), pt_num_(pts3d.shape(1)), frame_num_(pts2d.size())
    {
        assert(pts3d.shape(0) == 3);
        for(size_t i = 0;  i < frame_num_; i++)
        {
            assert(pts2d(i).shape(0) == 2);
            assert(pts2d(i).shape(1) == pts3d.shape(1));
            // pts2d_.push_back(Matrix<DualNumber<DType>>(pts2d(i)));
        }




        residual_ = Vector<DType>(frame_num_ * pt_num_ * 2);
        jacobian_pose_.resize(frame_num_);
        for(auto & jac : jacobian_pose_)
        {
            jac = Matrix<DType>({pt_num_ * 2, POSE_DOF});
        }
        jacobian_intrinsic_ = Matrix<DType>({frame_num_ * pt_num_ * 2, INTRINSIC_PARAM_NUM });
    }

    template<typename DTypeLocal=DType>
    Matrix<DTypeLocal> calcResidual(
        const Vector<RigidTransform<DTypeLocal, DIM>>& poses,
        const Camera<DTypeLocal>& cam) const
    {
        Vector<DTypeLocal> residual(frame_num_ * pt_num_ * 2);
        Camera<DTypeLocal> cam_shadow(cam);
        for(size_t frame_i = 0; frame_i < frame_num_; frame_i++)
        {
            cam_shadow.setPose(poses(frame_i));
            Matrix<DTypeLocal> pts2d = cam_shadow.project(pts3d_);
            Matrix<DTypeLocal> reprojection_error = pts2d - pts2d_(frame_i);
            residual.setBlock(frame_i * pt_num_ * 2, 0, reprojection_error(Row(0)).T());
            residual.setBlock(frame_i * pt_num_ * 2 + pt_num_, 0, reprojection_error(Row(1)).T());
        }

        return residual;
    }

    template<typename DTypeLocal=DType>
    Matrix<DTypeLocal> calcResidual(
        const Vector<RigidTransform<DTypeLocal, DIM>>& poses,
        const Camera<DTypeLocal>& cam,
        size_t frame_i) const
    {
        Vector<DTypeLocal> residual(pt_num_ * 2);
        Camera<DTypeLocal> cam_shadow(cam);

        cam_shadow.setPose(poses(frame_i));
        Matrix<DTypeLocal> pts2d = cam_shadow.project(pts3d_);
        Matrix<DTypeLocal> reprojection_error = pts2d - pts2d_(frame_i);
        residual.setBlock(0, 0, reprojection_error(Row(0)).T());
        residual.setBlock(pt_num_, 0, reprojection_error(Row(1)).T());

        return residual;
    }

    // RadialTangentialDistortion<DType>* pDistortion()
    // {
    //     if(!cam_.distortion()) return nullptr;
    //     return static_cast<RadialTangentialDistortion<DType>*>(cam_.distortion().get());
    // }

    // state_: (fx fy cx cy k1 k2 p1 p2) + (tx ty tx rx ry rz) * frame_num
    // update: vector t, rotate matrix r
    void update(const Vector<DType>& increment)
    {
        #if 1
        for(size_t frame_i = 0; frame_i < frame_num_; frame_i++)
        {
            Matrix<DType> r_update({3,1});
            for(size_t i = 0; i < DIM; i++)
            {
                poses_(frame_i).translation()(i) += increment(INTRINSIC_PARAM_NUM + frame_i * POSE_DOF + i);
                r_update(i, 0) = increment(INTRINSIC_PARAM_NUM + frame_i * POSE_DOF + i + DIM);
            }
            poses_(frame_i).rotation() =
            Rotation<DType, DIM>::fromMatrix(so::exp<DIM>(so::wedge(r_update))) *
            poses_(frame_i).rotation();
        }
        auto p_distortion = static_cast<RadialTangentialDistortion<DualNumber<DType>>*>(cam_.distortion().get());

        {
            Vector<DType> f = cam_.focalLength();
            Vector<DType> c = cam_.principalOffset();
            f(0) += increment(0);
            f(1) += increment(1);
            c(0) += increment(2);
            c(1) += increment(3);
            cam_.setFocalLength(f);
            cam_.setPrincipalOffset(c);
            p_distortion->k()(0) += increment(4);
            p_distortion->k()(1) += increment(5);
            p_distortion->p()(0) += increment(6);
            p_distortion->p()(1) += increment(7);
        }

        // update residual
        residual_ = calcResidual<DType>(poses_, cam_);

        // update jacobian
        Camera<DualNumber<DType>, 3> perturb_cam(cam_);
        auto p_perturb_distortion = static_cast<RadialTangentialDistortion<DualNumber<DType>>*>(perturb_cam.distortion().get());
        for(size_t frame_i = 0; frame_i < frame_num_; frame_i++)
        {
            // jacobian from translation
            for(size_t i = 0; i < DIM; i++)
            {
                poses_(frame_i).translation()(i)(1) = 1;
                Matrix<DualNumber<DType>> residual = calcResidual<DualNumber<DType>>(poses_, perturb_cam, frame_i);
                jacobian_pose_.at(frame_i)(Col(i)) = matrixAtPart(residual, 1);
                poses_(frame_i).translation()(i)(1) = 0;
            }

            // jacobian from rotation
            for(size_t i = 0; i < 3; i++)
            {
                Vector<DualNumber<DType>> rotvec{0,0,0};
                rotvec(i)(1) = 1;
                Rotation<DType, DIM> rotation_store = poses_(frame_i).rotation();
                poses_(frame_i).rotation() =
                    Rotation<DualNumber<DType>, DIM>::fromMatrix(so::exp<DIM>(so::wedge(rotvec)))
                    * poses_(frame_i).rotation();
                Matrix<DualNumber<DType>> residual = calcResidual<DualNumber<DType>>(poses_, perturb_cam, frame_i);
                jacobian_pose_.at(frame_i)(Col(i + DIM)) = matrixAtPart(residual, 1);
                poses_(frame_i).rotation() = rotation_store;
            }
        }


        // jacobian from distortion
        for(size_t i = 0; i < 2; i++)
        {
            Matrix<DualNumber<DType>> residual;
            Vector<DualNumber<DType>> f = cam_.focalLength();
            Vector<DualNumber<DType>> c = cam_.principalOffset();

            f(i)(1) = 1;
            perturb_cam.setFocalLength(f);
            residual = calcResidual<DualNumber<DType>>(poses_, perturb_cam);
            jacobian_intrinsic_(Col(i)) = matrixAtPart(residual, 1);
            perturb_cam.setFocalLength(cam_.focalLength());

            c(i)(1) = 1;
            perturb_cam.setPrincipalOffset(c);
            residual = calcResidual<DualNumber<DType>>(poses_, perturb_cam);
            jacobian_intrinsic_(Col(i + 2)) = matrixAtPart(residual, 1);
            perturb_cam.setPrincipalOffset(cam_.principalOffset());

            p_perturb_distortion->k()(i)(1) = 1;
            residual = calcResidual<DualNumber<DType>>(poses_, perturb_cam);
            jacobian_intrinsic_(Col(i + 4)) = matrixAtPart(residual, 1);
            p_perturb_distortion->k()(i)(1) = 0;

            p_perturb_distortion->p()(i)(1) = 1;
            residual = calcResidual<DualNumber<DType>>(poses_, perturb_cam);
            jacobian_intrinsic_(Col(i + 6)) = matrixAtPart(residual, 1);
            p_perturb_distortion->p()(i)(1) = 0;
        }
        #endif
    }

    void initialGuess(
        const Vector<RigidTransform<DType, DIM>>& poses_guess,
        const Camera<DType, DIM>& cam_guess)
    {
        // assert(poses_guess.size() == frame_num_);
        // for(size_t i = 0; i < frame_num_; i++)
        // {
        //     poses_.push_back(poses_guess(i));
        // }
        poses_ = poses_guess;
        cam_ = cam_guess;
        update(Vector<DType>::zeros( INTRINSIC_PARAM_NUM + poses_.size() * POSE_DOF ));
    }
    const Vector<RigidTransform<DualNumber<DType>, DIM>>& poses() const { return poses_; }
    // const RadialTangentialDistortion<DType>& distortion() const { return p_distortion_; }
    Camera<DType, DIM> camera() const { return cam_; }

    DType cost() const { return residual_.T().matmul(residual_)(0,0); }

    void solve(uint8_t step=20, uint8_t verbose=0)
    {
        Vector<DType> inc;
        for(size_t i = 0; i < step; i++)
        {
            if( cost() < 10 * eps<DType>()) break;
            inc = opt::solveSfMGaussNewton(jacobian_intrinsic_, jacobian_pose_, residual_);
            if(verbose > 0) std::cout << "\ni: " << i << std::endl;
            if(verbose > 0) std::cout << "cost: " << cost() << std::endl;

            // if(verbose > 2) ret += std::string("jac:\n") + mxm::to_string(problem.jac());
            // if(verbose > 3) ret += std::string("res:\n") + mxm::to_string(problem.res());
            update(inc);
        }
    }

private:
    Matrix<DualNumber<DType>> pts3d_;
    size_t pt_num_;
    size_t frame_num_;
    Camera<DualNumber<DType>, DIM> cam_;

    Vector<Matrix<DualNumber<DType>>> pts2d_;
    // std::shared_ptr<RadialTangentialDistortion<DType>> p_distortion_;
    Vector<RigidTransform<DualNumber<DType>, DIM>> poses_;
    Vector<DType> residual_;
    Matrix<DType> jacobian_intrinsic_;
    opt::BlockDiagMatrix<DType> jacobian_pose_;
};


} // namespace mxm



#endif // __CV_CALIBRATION_H__
