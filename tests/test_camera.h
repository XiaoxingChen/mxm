#if !defined(_TEST_CAMERA_H_)
#define _TEST_CAMERA_H_

#include "mxm/model_camera.h"

using namespace mxm;

inline void cameraTest1()
{
    Camera cam(RigidTrans::identity(3), Vec({100, 100}), Vec({2,2}));
    auto ray = cam.pixelRay({1,1});
    if((ray.direction() - Vec({0,0,1})).norm() > eps())
    {
        std::cout << ray.direction().T().str() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void cameraTest2()
{
    Camera cam(RigidTrans::identity(3), Vec({320,240}), Vec({640, 480}));
    Mat points({3,3}, {1, 1, 1, 0, 0, 1, -1,-1,1}, Mat::COL);

    Mat img = cam.project(points);
    Mat img_expected({2,3}, {640, 480, 320, 240, 0,0}, Mat::COL);

    if((img - img_expected).norm() > eps())
    {
        std::cout << "img: \n" << img.str() << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

inline void testCamera()
{
    cameraTest1();
    cameraTest2();
}


#endif // _TEST_CAMERA_H_
