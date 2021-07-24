#if !defined(_TEST_CAMERA_H_)
#define _TEST_CAMERA_H_
#include "test_config.h"

#if TEST_AVAILABLE_ALL
#include "mxm/model_camera.h"

using namespace mxm;

inline void cameraTest1()
{
    Camera<float> cam(RigidTrans::identity(3), Vec({100, 100}), Vec({2,2}));
    auto ray = cam.pixelRay({1,1});
    if((ray.direction() - Vec({0,0,1})).norm() > eps())
    {
        std::cout << mxm::to_string(ray.direction().T()) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void cameraTest2()
{
    Camera<float> cam(RigidTrans::identity(3), Vec({320,240}), Vec({640, 480}));
    Mat points({3,3}, {1, 1, 1, 0, 0, 1, -1,-1,1}, COL);

    Mat img = cam.project(points);
    Mat img_expected({2,3}, {640, 480, 320, 240, 0,0}, COL);

    if((img - img_expected).norm() > eps())
    {
        std::cout << "img: \n" << mxm::to_string(img) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testCamera()
{
    cameraTest1();
    cameraTest2();
    {
        Camera<double> cam;
        cam.setFocalLength({457.296, 458.654}).setResolution({480, 752}).setPrincipalOffset({248.375, 367.215});
        cam.setOrientation(Rotation<double>::fromAxisAngle({1,0,0}, M_PI_4));
        cam.setPosition({1,2,3});

        Matrix<size_t> px_src(fixRow(2), {100,100, 0,300, 234,123, 324,222, 335,123}, COL);

        auto px_pt = cam.pixelDirection(px_src) + cam.pose().translation();
        auto px_coord = cam.project(px_pt);
        double error(0);
        if(!isZero(px_src - px_coord, &error, 1e-12))
        {
            std::cout << "error: " << error << std::endl;
            std::cout << "eps: " << eps<double>() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        auto p_distor = Distortion<double>::radialTangential({-0.28340811, 0.07395907, 0.00019359, 1.76187114e-05});
        cam.setDistortion(p_distor);

        px_pt = cam.pixelDirection(px_src) + cam.pose().translation();
        px_coord = cam.project(px_pt);
        if(!isZero(px_src - px_coord, &error, 10))
        {
            std::cout << "Warning! Reprojection error: " << error << std::endl;
            // std::cout << mxm::to_string(px_src - px_coord) << std::endl;
            // throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }


    }
}
#else
inline void testCamera(){}
#endif

#endif // _TEST_CAMERA_H_
