#if !defined(_CV_BASIC_H_)
#define _CV_BASIC_H_

#include "cv_pixel.h"
#include "cv_kernel.h"
#include "interpolation.h"

namespace mxm
{
template<typename PType>
Matrix<PType>
resize(const Matrix<PType>& img, const Shape& shape, const std::string& strategy="bilinear")
{
    Matrix<PType> ret(shape);

    using DType = typename NormTraits<PType>::type;
    // std::cout << "k_w: " << k_w << ", k_h: " << k_h << std::endl;

    ret.traverse([&](auto i, auto j){
        if(strategy == std::string("bilinear"))
        {
            DType k_h = DType(img.shape(0) - 1) / (ret.shape(0) - 1);
            DType k_w = DType(img.shape(1) - 1) / (ret.shape(1) - 1);

            Matrix<PType> unit_mat({2,2});
            size_t i_0 = size_t(i * k_h);
            size_t j_0 = size_t(j * k_w);
            size_t i_1 = (i_0 == img.shape(0) - 1) ? i_0 : i_0 + 1;
            size_t j_1 = (j_0 == img.shape(1) - 1) ? j_0 : j_0 + 1;
            DType i_float = i * k_h - i_0;
            DType j_float = j * k_w - j_0;
            // std::cout << "i,j: " << i_int << "," << j_int << std::endl << std::flush;

            unit_mat(0,0) = img(i_0, j_0);
            unit_mat(0,1) = img(i_0, j_1);
            unit_mat(1,0) = img(i_1, j_0);
            unit_mat(1,1) = img(i_1, j_1);

            ret(i,j) = interp::bilinearUnitSquare(Vector<DType>({i_float, j_float}), unit_mat)(0,0);
            // todo, accelerate
        }else if(strategy == std::string("nearest"))
        {
            DType k_h = DType(img.shape(0)) / (ret.shape(0));
            DType k_w = DType(img.shape(1)) / (ret.shape(1));
            ret(i,j) = img(size_t(i * k_h + 0.), size_t(j * k_w + 0.));
        }else{
            throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
        }

        // ret(i,j) = interp::bilinearUnitSquare(Vector<DType>({0.5, 0.5}), unit_mat);
        // ret(i,j) = unit_mat(0,0);
    });
    return ret;
}

template<typename PType>
Matrix<PType>
resize(const Matrix<PType>& img, float rate, const std::string& strategy="bilinear")
{
    return resize(img, {size_t(img.shape(0) * rate), size_t(img.shape(1) * rate)}, strategy);
}

inline Matrix<size_t> nonMaximumSuppression(const Matrix<float>& heat_map, const Shape& block_shape=Shape({50,50}))
{
    // float thresh = 10;
    float thresh = mxm::sum(heat_map) / (heat_map.shape(0) * heat_map.shape(1));
    std::vector<size_t> local_maximums;
    for(size_t i = 0; i < heat_map.shape(0); i += block_shape[0])
    {
        for(size_t j = 0; j < heat_map.shape(1); j += block_shape[1])
        {
            auto maximal_idx = mxm::argMax(heat_map(Block(
                {i, std::min(i+block_shape[0], heat_map.shape(0))},
                {j, std::min(j+block_shape[1], heat_map.shape(1))})));

            if(heat_map(i+maximal_idx[0], j+maximal_idx[1]) < thresh) continue;

            local_maximums.push_back(i + maximal_idx[0]);
            local_maximums.push_back(j + maximal_idx[1]);

        }
    }
    return Matrix<size_t>(fixRow(2), std::move(local_maximums), Mat::COL);
}

inline Matrix<size_t> adaptiveNonMaximalSuppression(const Matrix<float>& heat_map)
{
    float thresh = mxm::sum(heat_map) / (heat_map.shape(0) * heat_map.shape(1));
    // thresh = 0;
    const size_t window_width = 5;
    const size_t half_width = window_width / 2;
    std::vector<size_t> coords;
    heat_map.traverse([&](auto i, auto j){
        if(i + window_width > heat_map.shape(0) || j + window_width > heat_map.shape(1)) return;
        if(heat_map(i+half_width, j+half_width) < thresh) return;
        if(mxm::argMax(heat_map(Block({i, i+window_width}, {j, j+window_width}))) == std::array<size_t, 2>{half_width, half_width})
        {
            coords.push_back(i+half_width);
            coords.push_back(j+half_width);
        }
    });
    // todo
    return Matrix<size_t>(fixRow(2), std::move(coords), Mat::COL);
}

inline Matrix<size_t> nmsGrid(const Matrix<float>& mat, const Shape& block_shape)
{
    std::vector<size_t> corners;
    for(size_t i = 0; i < mat.shape(0); i += block_shape[0])
    {
        for(size_t j = 0; j < mat.shape(1); j += block_shape[1])
        {
            corners.push_back(i);
            corners.push_back(j);
        }
    }
    return Matrix<size_t>(fixRow(2), std::move(corners), Mat::COL);
}

template<typename PType>
Matrix<PType> bilinear(const Matrix<float>& pts, const Matrix<PType>& img)
{
#if 0
    Matrix<PType> ret({pts.shape(1), 1});
    for(size_t i = 0;i < pts.shape(1); i++)
    {
        size_t x0(pts(0, i));
        size_t y0(pts(1, i));
        if(pts(0,i) + 1 >= img.shape(0) || pts(1,i) + 1 >= img.shape(1))
        {
            ret(i,0) = img(x0, y0);
            continue;
        }

        Matrix<PType> t = pts(Col(i)) - Matrix<float>({2,1}, {std::floor(pts(0, i)), std::floor(pts(1, i))});
        ret(i, 0) = interp::bilinearUnitSquare(t, img(Block({x0, x0+2},{y0, y0+2})));
    }
    return ret;
#else
    Matrix<PType> ret({pts.shape(1), 1});
    ret = interp::bilinearUnitSquare(pts, img);
    return ret;
#endif
}

template<typename PType> using ImagePyramid = std::vector<Matrix<PType>>;
template<typename PType>
std::vector<Matrix<PType>> gaussPyramid(const Matrix<PType>& img, size_t level, bool return_blured=true)
{
    using DType = typename NormTraits<PType>::type;
    std::vector<Matrix<PType>> pyramid;
    Matrix<PType> prev_img(img);
    for(size_t i = 0; i < level; i++)
    {
        auto blured = convolute(prev_img, kernel::gauss<DType>(3));
        prev_img = std::move(reduce(blured, kernel::average<DType>(2)));
        pyramid.push_back(std::move(blured));
    }
    return pyramid;
}

} // namespace mxm

#endif // _CV_BASIC_H_
