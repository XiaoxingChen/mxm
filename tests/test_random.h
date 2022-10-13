#if !defined(_TEST_RANDOM_H_)
#define _TEST_RANDOM_H_
#include "test_config.h"

#if TEST_AVAILABLE_RANDOM
#include "mxm/random.h"
#include "mxm/stats.h"
using namespace mxm;

// cov mat: diag{0.01, 0.04, 0.09}
// mean: {0,0,1}
Matrix<float> gaussianData01()
{
    Matrix<float> data({3, 100},
        {0.074678, 0.089794, -0.079578, -0.174893, 0.114849, 0.116078, 0.050221, -0.043998, 0.084360, 0.053618, 0.064195, 0.130798, -0.064995, -0.064926, -0.019389, -0.073063, -0.147800, -0.070540, -0.110595, 0.098409, 0.068942, 0.067242, -0.165465, -0.058896, 0.069746, -0.055833, -0.028598, -0.102073, -0.119206, 0.037265, -0.033150, 0.031803, -0.034244, -0.022900, 0.059597, -0.124129, 0.097870, 0.065079, 0.105512, 0.081861, 0.113410, 0.144086, 0.013386, -0.060453, -0.071908, -0.108322, -0.067986, 0.035513, 0.015830, 0.159110, 0.214225, 0.131257, 0.122079, -0.102439, -0.017687, -0.018661, 0.114281, -0.004656, 0.023987, 0.132209, -0.086739, 0.112669, 0.008727, -0.111474, 0.078892, -0.169443, -0.269868, 0.070571, -0.073534, -0.136648, -0.106449, -0.214329, -0.237239, -0.050980, -0.039303, -0.020569, -0.135250, -0.071593, 0.120580, 0.117533, 0.176392, 0.010753, 0.014947, 0.007419, 0.070004, -0.006251, -0.043432, -0.195877, 0.210981, 0.022418, 0.163364, -0.151808, -0.118890, 0.129483, 0.147191, 0.269179, 0.084686, -0.133734, 0.082787, 0.009141,
        0.089941, -0.173421, -0.299248, -0.304340, -0.236525, -0.000174, 0.225666, 0.056934, 0.106895, 0.021150, -0.229929, 0.036691, -0.009517, 0.050417, 0.040179, -0.221750, -0.015497, -0.334453, 0.111354, -0.008687, -0.051133, -0.063690, 0.209541, -0.292991, 0.069854, -0.172396, 0.242280, 0.097235, -0.030700, -0.328466, -0.253243, -0.094458, -0.219383, -0.045045, -0.023748, 0.118487, -0.039063, -0.124427, -0.117811, -0.079730, -0.046515, 0.294698, -0.129784, -0.409420, -0.233698, 0.124457, 0.123056, 0.054044, 0.118696, -0.082190, 0.136687, -0.659527, 0.243109, 0.076856, -0.047415, 0.321302, 0.027839, 0.094665, 0.253163, -0.047303, 0.224743, 0.210043, -0.091875, 0.153181, 0.228630, 0.592785, -0.038725, -0.178630, 0.201075, 0.040964, 0.333992, -0.087168, 0.079767, 0.136593, 0.127142, 0.185122, 0.013088, 0.202200, -0.147079, -0.012407, -0.019080, 0.368287, 0.246718, 0.227243, -0.224230, 0.027313, 0.050938, -0.005168, -0.055115, -0.290466, 0.031644, -0.397755, -0.031716, 0.200747, 0.276541, -0.029280, -0.243402, -0.113771, 0.272885, 0.306044,
        0.922272, 1.480478, 0.550312, 1.052430, 1.035779, 0.909393, 1.137454, 1.056695, 1.118492, 0.289048, 1.013259, 0.800652, 0.809123, 0.948626, 1.519379, 0.770544, 1.210945, 1.043605, 1.028101, 1.361236, 1.078306, 0.792412, 1.504428, 0.789392, 1.052560, 1.082252, 0.741081, 1.355869, 1.394962, 1.204310, 1.282769, 1.248510, 0.823084, 1.062731, 1.258137, 1.088326, 0.859938, 0.499365, 0.960477, 1.472987, 0.762486, 0.723305, 0.685811, 0.867411, 1.286611, 0.990282, 1.003084, 1.275928, 0.398766, 1.054127, 0.998780, 1.098587, 0.873495, 0.927715, 0.274478, 1.004428, 1.131307, 0.929267, 0.290341, 1.593756, 1.230139, 1.222738, 0.970087, 1.026175, 1.099059, 1.281591, 1.270393, 1.283363, 0.772539, 0.686656, 1.058493, 0.343471, 1.335673, 1.088399, 1.343335, 0.911611, 0.825859, 0.525193, 0.404265, 1.692690, 0.745036, 0.831340, 1.071557, 0.654193, 0.956754, 1.980213, 1.289231, 1.040404, 1.317409, 1.115859, 0.868219, 0.843912, 0.794532, 1.119370, 0.924340, 0.862714, 0.834807, 0.711095, 1.094230, 1.239817});

    return data;
}

void testStatistic02()
{
    Matrix<float> data(fixRow(3), {
        0.78472045, 0.44620021, 0.15966676, 0.63083727, 0.454521, 0.67257604, 0.94173523, 0.12610781, 0.12999695, 0.88745649,
       0.61714868, 0.76269721, 0.40183527, 0.64904284, 0.52523561, 0.50782472, 0.89589038, 0.86108194, 0.12262409, 0.37116328,
       0.18776345, 0.05760265, 0.89946426, 0.63555678, 0.1758316, 0.39483145, 0.89101437, 0.62559045, 0.50012674, 0.4658263 });

    Vector<float> expected_mean({0.52338182, 0.5714544 , 0.4833608});
    Matrix<float> expected_cov({3,3}, {
        0.09613769,  0.02057896, -0.00781674,
        0.02057896,  0.05670503,  0.00395378,
       -0.00781674,  0.00395378,  0.08381827});

    auto mean = findMean(data);
    auto cov = findCovariance(data, mean);

    float error;
    if(  !isZero(expected_cov - cov, &error, 0.0001))
    {

        std::cout << "expected_cov: \n" << mxm::to_string(expected_cov) << std::endl;
        std::cout << "cov: \n" << mxm::to_string(cov) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(  !isZero(mean - expected_mean, &error, 0.0001))
    {
        std::cout << "mean: \n" << mxm::to_string(mean) << std::endl;
        std::cout << "expected_mean: \n" << mxm::to_string(expected_mean) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }


}

void testGaussian()
{
    Vector<float> mean({1,2,3});
    auto real_cov = mxm::diagonalMatrix(Vector<float>{0.0001, 0.01, 0.04});
    auto data = random::gaussian(mean, real_cov, 100);
    auto generated_cov = findCovariance(data, mean);
    auto generated_mean = findMean(data);

    float error;
    if(  !isZero(real_cov - generated_cov, &error, 0.01))
    {

        std::cout << "real_cov: \n" << mxm::to_string(real_cov) << std::endl;
        std::cout << "generated_cov: \n" << mxm::to_string(generated_cov) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(  !isZero(mean - generated_mean, &error, 0.05))
    {
        std::cout << "mean: \n" << mxm::to_string(mean) << std::endl;
        std::cout << "generated_mean: \n" << mxm::to_string(generated_mean) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
}

void testStatistic01()
{
    Vector<float> mean({0,0,1});
    auto real_cov = mxm::diagonalMatrix(Vector<float>{0.01, 0.04, 0.09});

    auto data = gaussianData01();
    auto generated_cov = findCovariance(data, mean);
    auto generated_mean = findMean(data);

    // std::cout << "data: \n" << mxm::to_string(data) << std::endl;

    float error;
    if(  !isZero(real_cov - generated_cov, &error, 0.01))
    {

        std::cout << "real_cov: \n" << mxm::to_string(real_cov) << std::endl;
        std::cout << "generated_cov: \n" << mxm::to_string(generated_cov) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

    if(  !isZero(mean - generated_mean, &error, 0.01))
    {
        std::cout << "mean: \n" << mxm::to_string(mean) << std::endl;
        std::cout << "generated_mean: \n" << mxm::to_string(generated_mean) << std::endl;
        std::cout << "error: " << error << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }
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
    testStatistic01();
    testStatistic02();
}
#else
void testRandom(){}
#endif
#endif // _TEST_RANDOM_H_
