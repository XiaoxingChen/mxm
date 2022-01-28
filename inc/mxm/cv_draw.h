#if !defined(_CV_DRAW_H_)
#define _CV_DRAW_H_

#include "cv_basic.h"
#include <iostream>

namespace mxm
{

Matrix<size_t> pixelsOnLine(const Matrix<float>& pts, float width=5)
{
    Matrix<float> delta = pts(Col(1)) - pts(Col(0));
    Matrix<float> dir_inc(delta.shape(), {0.f + delta(0,0), 0.f + delta(1,0)});

    size_t main_ax = abs(delta(0,0)) > abs(delta(1,0)) ? 0 : 1;
    size_t sub_ax = (main_ax == 0) ? 1 : 0;

    dir_inc *= (1.f / abs(delta(main_ax,0)));
    float sub_width = abs(dir_inc.normalized()(main_ax, 0)) * width;

    Matrix<float> center({2,1}, {0.f + pts(0,0), 0.f + pts(1,0)});
    std::vector<size_t> ret_data;
    // std::cout << "main ax: " << main_ax << std::endl;
    // std::cout << "dir_inc: " << mxm::to_string(dir_inc.T()) << std::endl;

    for(size_t t = 0; t < abs(delta(main_ax,0)); t++)
    {
        for(float u = center(sub_ax, 0) - sub_width; u < center(sub_ax, 0) + sub_width; u += 1)
        {
            if(main_ax == 0)
            {
                ret_data.push_back(size_t(center(main_ax, 0) + 0.5));
                ret_data.push_back(size_t(u + 0.5));
            }else
            {
                ret_data.push_back(size_t(u + 0.5));
                ret_data.push_back(size_t(center(main_ax, 0) + 0.5));
            }
        }
        center += dir_inc;
        // std::cout << "t: " << t << std::endl;
    }

    Matrix<size_t> ret({2, ret_data.size() / 2}, std::move(ret_data), COL);
    return ret;
}

template<typename PType>
void scatter(Matrix<PType>& p, const Matrix<float>& pts, const PType& color=PType::black(), float width=5)
{
    // auto color = PType({0,0,128});
    for(size_t i = 0; i < pts.shape(1); i++)
    {
        if(pts(0, i) < 0 || pts(0, i) > p.shape(0) || pts(1, i) < 0 || pts(1, 0) > p.shape(1)) continue;

        size_t start_x = std::max<float>(pts(0, i)-width*0.5, 0) + 0.5f;
        size_t start_y = std::max<float>(pts(1, i)-width*0.5, 0) + 0.5f;

        if(start_x > p.shape(0) || start_y > p.shape(1)) continue;

        size_t width_x = std::min<float>(width, p.shape(0) - 1 - start_x) + 0.5f;
        size_t width_y = std::min<float>(width, p.shape(1) - 1 - start_y) + 0.5f;
        p.setBlock(start_x, start_y, Matrix<PType>::zeros({width_x, width_y}) + color);
    }
}

template<typename PType>
void
plot(Matrix<PType>& p, const Matrix<float>& pts, const PType& color=PType::black(), float width=2)
{
    // auto color = PType({0,0,0.8});
    #if 0
    for(size_t i = 0; i < pts.shape(1); i++)
    {
        size_t start_x = std::max<float>(pts(0, i)-width*0.5, 0) + 0.5f;
        size_t start_y = std::max<float>(pts(1, i)-width*0.5, 0) + 0.5f;

        if(start_x > p.shape(0) || start_y > p.shape(1)) continue;

        size_t width_x = std::min<float>(width, p.shape(0) - 1 - start_x) + 0.5f;
        size_t width_y = std::min<float>(width, p.shape(1) - 1 - start_y) + 0.5f;
        p.setBlock(start_y, start_x, Matrix<PType>::zeros({width_y, width_x}) + color);
    }
    #endif
    #if 1

    for(size_t i = 0; i < pts.shape(1) - 1; i++)
    {
        auto pxs = pixelsOnLine(pts(Block({0, end()}, {i, i+2})), width);
        for(size_t j = 0; j < pxs.shape(1); j++)
        {
            if(pxs(0,j) >= p.shape(0)) continue;
            if(pxs(1,j) >= p.shape(1)) continue;
            p(pxs(0,j), pxs(1,j)) = color;
        }
    }
    #endif

}

template<typename PType> void
drawLines(Matrix<PType>& p, const Matrix<float>& vertex, const Matrix<size_t>& indices, const PType& color=PType::black(), float width=2)
{
    for(size_t i = 0; i < indices.shape(1); i++)
    {
        auto idx0 = indices(0, i);
        auto idx1 = indices(1, i);
        auto pts = Matrix<float>({2,2}, {vertex(0, idx0), vertex(1, idx0), vertex(0, idx1), vertex(1, idx1)}, COL);
        pts(0, 0) = std::max(std::min(float(p.shape(0)), pts(0, 0)), 0.f);
        pts(0, 1) = std::max(std::min(float(p.shape(0)), pts(0, 1)), 0.f);

        pts(1, 0) = std::max(std::min(float(p.shape(1)), pts(1, 0)), 0.f);
        pts(1, 1) = std::max(std::min(float(p.shape(1)), pts(1, 1)), 0.f);

        auto pxs = pixelsOnLine(pts, width);
        for(size_t j = 0; j < pxs.shape(1); j++)
        {
            if(pxs(0,j) >= p.shape(0)) continue;
            if(pxs(1,j) >= p.shape(1)) continue;
            p(pxs(0,j), pxs(1,j)) = color;
        }
    }
}

} // namespace mxm



#endif // _CV_DRAW_H_
