#if !defined(_CV_OPTICAL_FLOW_H_)
#define _CV_OPTICAL_FLOW_H_

#include "cv_basic.h"
#include "linalg_utils.h"
#include "linalg_solve.h"
#include "cv_kernel.h"

namespace mxm
{

// input:
// x,y: center of the window
// return: all pixel coordinates in the window
template<typename DType=size_t>
inline Matrix<DType> windowGrids(DType x, DType y, const Shape& shape_limit, size_t window_width)
{
    std::vector<DType> ret_data;
    size_t half_w = window_width / 2;
    for(DType i = x > half_w ? x - half_w : 0 ; i < std::min(x + half_w + 1, DType(shape_limit[0])); i++)
    {
        for(DType j = y > half_w ? y - half_w : 0; j < std::min(y + half_w + 1, DType(shape_limit[1])); j++)
        {
            ret_data.push_back(i);
            ret_data.push_back(j);
        }
    }
    return Matrix<DType>(fixRow(2), std::move(ret_data), Mat::COL);
}

// return: updated point position
// negative position if point lost tracking.
inline Matrix<float> lkOpticalFlow(
    const Matrix<float>& img_curr,
    const Matrix<float>& img_next,
    const Matrix<float>& pts_curr,
    const Matrix<float>& pts_init,
    size_t window_width=5)
{
    Matrix<float> sobel_x = convolute(img_curr, kernel::sobel());
    Matrix<float> sobel_y = convolute(img_curr, kernel::sobel().T());

    const size_t max_it = 10;
    const float tol = 1e-2;
    const bool verbose = 0;
    const float max_move = 4;
    const Matrix<float> coord_lost = Matrix<float>::ones({2,1}) * -1;

    Matrix<float> pts_next(pts_curr.shape());

    for(size_t pt_idx = 0; pt_idx < pts_curr.shape(1); pt_idx++)
    {
        if(pts_init(0, pt_idx) < 0)
        {
            pts_next(Col(pt_idx)) = coord_lost;
            continue;
        }

        Matrix<float> grids_curr = windowGrids<float>(pts_curr(0, pt_idx), pts_curr(1, pt_idx), img_curr.shape(), window_width);
        Matrix<float> pt_next = pts_init(Col(pt_idx));
        float prev_inc = 10.;
        if(grids_curr.shape(1) != window_width * window_width)
        {
            pts_next(Col(pt_idx)) = coord_lost;
            continue;
        }
        bool track_lost = false;

        for(size_t _i = 0; _i < max_it; _i++)
        {
            Matrix<float> grids_next = windowGrids<float>(pt_next(0,0), pt_next(1,0), img_next.shape(), window_width);
            if(grids_next.shape(1) != window_width * window_width) {track_lost = true; break;}
            Matrix<float> mat_a({grids_curr.shape(1), 2});
            Matrix<float> vec_b({grids_next.shape(1), 1});

            mat_a(Col(0)) = bilinear(grids_curr, sobel_x);
            mat_a(Col(1)) = bilinear(grids_curr, sobel_y);
            vec_b(Col(0)) = bilinear(grids_curr, img_curr) - bilinear(grids_next, img_next);

            Matrix<float> solve_a(mat_a.T().matmul(mat_a));
            Matrix<float> solve_b(mat_a.T().matmul(vec_b));

            Matrix<float> increment = qr::solve(solve_a, solve_b);
            if(increment.norm() > prev_inc) { std::cout << "inc inc" << std::endl; break; }
            pt_next += increment;
            if(verbose) std::cout << "increment: " << increment.norm() << std::endl;
            if(increment.norm() < tol) break;
        }
        if(verbose) std::cout << "done" << std::endl;
        if(!track_lost && (pt_next - pts_init(Col(pt_idx))).norm() > max_move)
        {
            if(verbose) std::cout << "pt jump" << std::endl;
            track_lost = true;
        }
        pts_next(Col(pt_idx)) = track_lost ? coord_lost : pt_next;
    }
    return pts_next;
}

inline Matrix<float> lkOpticalFlow(
    const Matrix<float>& img_curr,
    const Matrix<float>& img_next,
    const Matrix<float>& pts_curr,
    size_t window_width=5)
{
    return lkOpticalFlow(img_curr, img_next, pts_curr, pts_curr, window_width);
}

inline Matrix<float> pyramidLKOpticalFlow(
    const std::vector<Matrix<float>>& img_curr,
    const std::vector<Matrix<float>>& img_next,
    const Matrix<float>& pts_curr,
    size_t window_width=5)
{
    Matrix<float> pts_next = pts_curr * (1./ float(1 << (img_curr.size() - 1)));
    for(int layer = img_curr.size() - 1; layer > 0; layer--)
    {
        float scale = (1. / float(1 << layer));
        pts_next = lkOpticalFlow(img_curr.at(layer), img_next.at(layer), pts_curr * scale, pts_next, window_width);
        pts_next *= 2;
    }
    return pts_next;
}
} // namespace mxm


#endif // _CV_OPTICAL_FLOW_H_
