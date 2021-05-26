#if !defined(_TEST_PIXEL_H_)
#define _TEST_PIXEL_H_
#include "test_config.h"

#if TEST_AVAILABLE_ALL
#include <iostream>
#include "mxm/cv_basic.h"
#include "mxm/linalg.h"
#include "mxm/cv_kernel.h"


using namespace mxm;

inline void testPixel()
{
    Pixel color({1,1,1});
    if(color.rU8() != 0xff)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    Pixel px2(color);
    if(px2.rU8() != 0xff)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    std::vector<Pixel> img(10, px2);
    if(img.back().rU8() != 0xff)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    std::vector<Pixel> img2(img);
    img.insert(img.end(), img2.begin(), img2.begin() + 5);
    if(img.back().rU8() != 0xff)
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));


    {
        auto mat_a = random::uniform<Pixel>({3,3});
        if(mat_a(0,0).r() > 1)
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }

}

inline void testImageResize()
{
    {// case 1
        Matrix<Pixel> img({2,2});
        img(0,0) = Pixel::black();
        img(0,1) = Pixel::white();
        img(1,0) = Pixel::black();
        img(1,1) = Pixel::white();

        auto img3x3 = resize(img, Shape({3,3}));
        if(img3x3(1,1) != Pixel({.5, .5, .5}))
        {
            std::cout << to_string(img3x3,2) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }

    {//case 2
        Matrix<float> img({2,2},{0,1,2,3});
        Matrix<float> expected = Matrix<float>::ones({8,8});
        expected(Block({0,4},{0,4}))*=0;
        expected(Block({0,4},{4,8}))*=1;
        expected(Block({4,8},{0,4}))*=2;
        expected(Block({4,8},{4,8}))*=3;

        if(mxm::resize(img, 4, "nearest") != expected)
        {
            std::cout << to_string(mxm::resize(img, 4, "nearest"),2) << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }
    }
}

inline void testPixelMemory()
{
    std::vector<Pixel> pxs;
    pxs.push_back(Pixel({0.5,0.25,0.125}));
    pxs.push_back(Pixel({0.5,0.25,0.125}));
    pxs.push_back(Pixel({0.5,0.25,0.125}));

    for(size_t i = 0;i < 9; i+= Pixel::size())
    {
        if(
            ((float*)pxs.data())[i+0] != 0.5 ||
            ((float*)pxs.data())[i+1] != 0.25 ||
            ((float*)pxs.data())[i+2] != 0.125)
        {
            std::cout << ((float*)pxs.data())[i+0] << ","
                << ((float*)pxs.data())[i+1] << ","
                << ((float*)pxs.data())[i+2] << std::endl;
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

    }
}

inline void testKernels()
{
    if(0){//case 01
        std::cout << mxm::to_string(kernel::gauss<float>(3)) << std::endl;
        std::cout << mxm::sum(kernel::gauss<float>(3)) << std::endl;
    }
}

inline void testCvBasic()
{
    testPixel();
    testImageResize();
    testPixelMemory();
    testKernels();
}
#else
inline void testCvBasic(){}
#endif
#endif // _TEST_PIXEL_H_
