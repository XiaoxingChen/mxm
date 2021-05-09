#if !defined(_RANDOM_H_)
#define _RANDOM_H_

#include <functional>
#include <random>
#include "linalg_mat.h"
#include "cv_pixel.h"

namespace mxm
{
namespace random
{

template<typename DType>
typename std::enable_if<std::is_floating_point<DType>::value, DType>::type
uniform()
{
    static std::uniform_real_distribution<DType> distribution(0.0, 1.0);
    static std::mt19937 generator;
    static std::function<DType()> rand_generator = std::bind(distribution, generator);
    return rand_generator();
}

template<typename DType>
Matrix<DType> uniform(const Shape& shape)
{
    Matrix<DType> ret(shape);
    ret.traverse([&](auto i, auto j) { ret(i,j) = uniform<DType>(); });
    return ret;
}

template<typename DType>
inline Matrix<DType> unitSphere(size_t dim=3, size_t num=1)
{
    Matrix<DType> p({dim, num});
    for(size_t i = 0; i < num; i++)
    {
        do{
            p(Col(i)) = uniform<DType>({dim, 1});
        } while(p(Col(i)).norm() > 1);
    }

    return p;
}

} // namespace ramdom

} // namespace mxm



#endif // _RANDOM_H_
