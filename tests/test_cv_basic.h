#if !defined(_TEST_PIXEL_H_)
#define _TEST_PIXEL_H_

#include <iostream>
#include "mxm/cv_basic.h"
#include "mxm/linalg.h"


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

    // img3x3.traverse([&](auto i, auto j){
    //     std::cout << to_string(img3x3(i,j)) << std::endl;
    // });
    // std::cout << "TODO: to_string(const Matrix<Pixel>&) doesn't work well!" << std::endl;
    // std::cout << __FILE__ << ":" << __LINE__ << std::endl;
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

inline void testCvBasic()
{
    testPixel();
    testImageResize();
    testPixelMemory();
}

#endif // _TEST_PIXEL_H_
