#if !defined(_RANDOM_H_)
#define _RANDOM_H_

#include <functional>
#include <random>
#include "linalg_mat.h"
#include "random_forward_declaration.h"

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

template<typename DType>
Matrix<DType> gaussian(const Vector<DType>& mean, const Matrix<DType>& cov, size_t num)
{
    assert(mean.size() == cov.shape(0) && cov.square());
    std::default_random_engine generator;

    Matrix<DType> ret({mean.size(), num});
    auto val_vec = symmetricEig(cov);
    for(size_t i = 0; i < mean.size(); i++)
    {
        std::normal_distribution<DType> distribution(mean(i), sqrt(val_vec[0](i,0)));
        for(size_t j = 0; j < num; j++)
        {
            ret(i,j) = distribution(generator);
            // ret(i,j) = 0;
        }
    }

    return val_vec[1].matmul(ret);
    // return (ret);
}

template<typename DType>
size_t weightedSample(const Vector<DType>& weights)
{
    if(weights.size() < 2) return 0;

    Vector<DType> accumu_weights = weights;
    size_t n = weights.size();
    for(size_t i = 1; i < n; i++)
    {
        accumu_weights(i) += accumu_weights(i-1);
    }
    if(accumu_weights(n - 1) > DType(1) + eps<DType>())
    {
        return weightedSample(Vector<DType>(normalized(weights)));
    }
    DType rand_val = uniform<DType>();
    for(size_t i = 0; i < n; i++)
    {
        if(rand_val < accumu_weights(i)) return i;
    }
    return n-1;
}

} // namespace ramdom

} // namespace mxm



#endif // _RANDOM_H_
