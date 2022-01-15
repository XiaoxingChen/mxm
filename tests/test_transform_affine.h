#if !defined(_TEST_TRANSFORM_AFFINE_H_)
#define _TEST_TRANSFORM_AFFINE_H_

#include "mxm/transform_affine.h"
#include "mxm/rotation.h"
#include "mxm/geometry_ray.h"
using namespace mxm;

inline AffineTransform<float, 3> affineTransformTestData0()
{
    AffineTransform<float, 3> tf({1,1,1}, Rotation<float, 3>::fromAxisAngle({1,0,0}, 0.5), {1,2,3});
    return tf;
}

inline AffineTransform<float, 3> affineTransformTestData1()
{
    Matrix<float> mat({4,4},{1,2,3,10, 3,2,1,10, 1,3,9,10, 0,0,0,1}, ROW);
    auto tf = AffineTransform<float, 3>::fromMatrix(mat);
    return tf;
}

inline void testAffineTransform()
{
    {
        const size_t DIM = 3;
        AffineTransform<float, 3> tf({1,1,1}, Rotation<float, DIM>::fromAxisAngle({1,0,0}, 0.5), {1,2,3});

        float error;
        if(!isIdentity((tf.inv() * tf).asMatrix(), &error, 1e-5))
        {
            std::cout << mxm::to_string((tf.inv() * tf).asMatrix()) << std::endl;
            auto mat = tf.asMatrix();
            std::cout << mxm::to_string(mxm::inv(mat).matmul(mat)) << std::endl;

            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        tf = affineTransformTestData1();

        if(!isIdentity((tf.inv() * tf).asMatrix(), &error, 1e-5))
        {
            std::cout << mxm::to_string((tf.inv() * tf).asMatrix()) << std::endl;
            auto mat = tf.asMatrix();

            std::cout << mxm::to_string(mxm::inv(mat).matmul(mat)) << std::endl;
            std::cout << mxm::to_string(mat) << std::endl;
            std::cout << mxm::to_string(mxm::inv(mat)) << std::endl;

            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {
        const size_t DIM = 3;

        auto tf = affineTransformTestData1();
        if(!isZero(tf.asMatrix() - AffineTransform<float, DIM>::fromMatrix(tf.asMatrix()).asMatrix(), nullptr, 1e-5))
        {
            std::cout << mxm::to_string(tf.asMatrix()) << std::endl;

            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

    {
        auto tf = affineTransformTestData1();
        Ray ray({1,2,3}, {0,1,2});

        auto inv_ray = apply(tf, ray);
        auto expect_ray = apply(tf.inv(), inv_ray);

        if(! isZero(ray.origin() - expect_ray.origin(), nullptr, 1e-5))
        {
            std::cout << mxm::to_string(expect_ray.origin()) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(! isZero(ray.direction() - expect_ray.direction(), nullptr, 1e-5))
        {
            std::cout << mxm::to_string(expect_ray.direction()) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(abs(ray.tMin() - expect_ray.tMin()) > 1e-5)
        {
            std::cout << expect_ray.tMin() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        if(abs(ray.tMax() - expect_ray.tMax()) > 1e-2)
        {
            std::cout << ray.tMax() - expect_ray.tMax() << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }


    }



}

#endif // _TEST_TRANSFORM_AFFINE_H_
