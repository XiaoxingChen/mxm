#include "test_config.h"
#include "mxm/linalg_dual_number.h"

using namespace mxm;

template<typename DType>
DType testFunction01(DType x)
{
    return 3.f*x*x + 2.f*x;
}

void dualNumberTest01()
{
    float x = 2.f;
    float dfdx = testFunction01(DualNumber<float>{x, 1})(1);
    if(abs(dfdx - 14) > eps())
    {
        std::cout << "dfdx: " << dfdx << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

void testDualNumber()
{
    dualNumberTest01();
}
