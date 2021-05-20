#if !defined(_CV_CORNER_H_)
#define _CV_CORNER_H_

#include "cv_pixel.h"
#include "linalg_utils.h"

namespace mxm
{
Matrix<float> harrisCornernessMap(const Matrix<float>& src, size_t window_width=5, float k=0.06)
{
    Matrix<float> sobel = Matrix<float>({3,3}, {-1,0,1, -2,0,2, -1,0,1});
    auto sobel_x = convolute(src, sobel.T());
    auto sobel_y = convolute(src, sobel);

    Matrix<float> result(src.shape());
    result.traverse([&](auto i, auto j){
        if(i+window_width >= sobel_x.shape(0) || j+window_width >= sobel_x.shape(1)) return;

    #if 0
        Matrix<float> M({2,2});
        auto Ix = sobel_x(Block({i, i+window_width}, {j, j+window_width}));
        auto Iy = sobel_y(Block({i, i+window_width}, {j, j+window_width}));
        M(0,0) = mxm::sum(Ix * Ix);
        M(1,1) = mxm::sum(Iy * Iy);
        M(0,1) = mxm::sum(Ix * Iy);
    #else
        Matrix<float> M = Matrix<float>::zeros({2,2});
        for(size_t u = i; u < i+window_width; u++)
        {
            for(size_t v = j; v < j+window_width; v++)
            {
                float ix = sobel_x(u,v);
                float iy = sobel_y(u,v);
                M(0,0) += ix*ix;
                M(1,1) += iy*iy;
                M(0,1) += ix*iy;
            }
        }
    #endif
        M(1,0) = M(0,1);

        float det = M.det();
        float tr = M.trace();

        result(i+window_width/2, j+window_width/2) = det - k * tr*tr;
    });
    return result;
}
} // namespace mxm

#endif // _CV_CORNER_H_
