#if !defined(_TEST_TRANSFORM_AFFINE_H_)
#define _TEST_TRANSFORM_AFFINE_H_

#include "mxm/transform_affine.h"
#include "mxm/rotation.h"
using namespace mxm;

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
    }

    {
        const size_t DIM = 3;
        Matrix<float> mat({4,4},{1,2,3,11, 4,5,6,12, 7,8,9,13, 0,0,0,1}, ROW);
        auto tf = AffineTransform<float, DIM>::fromMatrix(mat);
        if(!isZero(tf.asMatrix() - mat, nullptr, 1e-5))
        {
            std::cout << mxm::to_string(tf.asMatrix()) << std::endl;

            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }

}

#endif // _TEST_TRANSFORM_AFFINE_H_
