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

}

inline void testImageResize()
{
    Matrix<Pixel> img({2,2});
    img(0,0) = Pixel::black();
    img(0,1) = Pixel::white();
    img(1,0) = Pixel::black();
    img(1,1) = Pixel::white();

    auto img3x3 = resize(img, Shape({3,3}));

    // img3x3.traverse([&](auto i, auto j){
    //     std::cout << to_string(img3x3(i,j)) << std::endl;
    // });
    std::cout << "TODO: to_string(const Matrix<Pixel>&) doesn't work well!" << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
}

inline void testCvBasic()
{
    testPixel();
    testImageResize();
}

#endif // _TEST_PIXEL_H_
