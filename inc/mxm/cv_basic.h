#if !defined(_CV_BASIC_H_)
#define _CV_BASIC_H_

#include "cv_pixel.h"
#include "interpolation.h"

namespace mxm
{
template<typename DType, size_t NChannel>
Matrix<PixelType<DType, NChannel>>
resize(const Matrix<PixelType<DType, NChannel>>& img, const Shape& shape)
{
    Matrix<PixelType<DType, NChannel>> ret(shape);

    FloatType k_h = FloatType(img.shape(0) - 1) / (ret.shape(0) - 1);
    FloatType k_w = FloatType(img.shape(1) - 1) / (ret.shape(1) - 1);
    std::cout << "k_w: " << k_w << ", k_h: " << k_h << std::endl;

    ret.traverse([&](auto i, auto j){
        Matrix<PixelType<DType, NChannel>> unit_mat({2,2});
        size_t i_0 = size_t(i * k_h);
        size_t j_0 = size_t(j * k_w);
        size_t i_1 = (i_0 == img.shape(0) - 1) ? i_0 : i_0 + 1;
        size_t j_1 = (j_0 == img.shape(1) - 1) ? j_0 : j_0 + 1;
        FloatType i_float = i * k_h - i_0;
        FloatType j_float = j * k_w - j_0;
        // std::cout << "i,j: " << i_int << "," << j_int << std::endl << std::flush;

        unit_mat(0,0) = img(i_0, j_0);
        unit_mat(0,1) = img(i_0, j_1);
        unit_mat(1,0) = img(i_1, j_0);
        unit_mat(1,1) = img(i_1, j_1);

        ret(i,j) = interp::bilinearUnitSquare(Vector<DType>({i_float, j_float}), unit_mat);
        // ret(i,j) = interp::bilinearUnitSquare(Vector<DType>({0.5, 0.5}), unit_mat);
        // ret(i,j) = unit_mat(0,0);
    });
    return ret;
}

} // namespace mxm

#endif // _CV_BASIC_H_
