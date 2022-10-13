#if !defined(__TEST_DSP_FFT__)
#define __TEST_DSP_FFT__

#include "mxm/dsp_fft.h"


namespace mxm{

inline void testDspFFT()
{
    {
        auto mat4 = mxm::dftMatrix(4);
        Matrix<Complex<float>> expected({4,4}, {
            1.f,1,1,1,
            1,{0,-1},-1,{0,1},
            1,-1,1,-1,
            1,{0,1},-1,{0,-1}
        }, ROW);
        expected *= 0.5;
        if(!isZero(mat4 - expected, nullptr, 1e-6))
        {
            std::cout << "mat4: " << mxm::to_string(mat4) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}
}
#endif // __TEST_DSP_FFT__
