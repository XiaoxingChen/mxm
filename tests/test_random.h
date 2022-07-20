#if !defined(_TEST_RANDOM_H_)
#define _TEST_RANDOM_H_
#include "test_config.h"

#if TEST_AVAILABLE_RANDOM
#include "mxm/random.h"
using namespace mxm;

void testGaussian()
{
    Vector<float> residual({0.1, 0.2, 0.3});
    Vector<float> mean({0,0,0});
    auto real_cov = residual.matmul(residual.T());
    auto data = random::gaussian(mean, real_cov, 100);
    auto cov = findCovariance(data, mean);
    // std::cout << "real_cov: " << mxm::to_string(real_cov) << std::endl;
    // std::cout << "cov: " << mxm::to_string(cov) << std::endl;
}

void testRandom()
{
    for(size_t dim = 2; dim < 5; dim++)
    {
        auto vec = random::unitSphere<FloatType>(dim);
        if(vec.shape(0) != dim)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    testGaussian();
}
#else
void testRandom(){}
#endif
#endif // _TEST_RANDOM_H_
