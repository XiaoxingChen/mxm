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

inline Matrix<float> lkOpticalFlow(
    const Matrix<float>& img_curr,
    const Matrix<float>& img_next,
    const Matrix<float>& pts_curr,
    size_t window_width=5)
{
    Matrix<float> sobel_x = convolute(img_curr, kernel::sobel());
    Matrix<float> sobel_y = convolute(img_curr, kernel::sobel().T());

    const size_t max_it = 40;
    float tol = 1e-3;

    Matrix<float> pts_next(pts_curr.shape());

    for(size_t pt_idx = 0; pt_idx < pts_curr.shape(1); pt_idx++)
    {
        Matrix<size_t> grids_curr = windowGrids<size_t>(pts_curr(0, pt_idx) + 0.5, pts_curr(1, pt_idx) + 0.5, img_curr.shape(), window_width);
        Matrix<float> pt_next = pts_curr(Col(pt_idx));
        float prev_inc = 10.;
        // if(grids_curr.shape(1) != window_width * window_width) continue;
        for(size_t _i = 0; _i < max_it; _i++)
        {
            Matrix<float> grids_next = windowGrids<float>(pt_next(0,0), pt_next(1,0), img_next.shape(), window_width);
            Matrix<float> mat_a({grids_curr.shape(1), 2});
            Matrix<float> vec_b({grids_next.shape(1), 1});
            for(size_t grid_idx = 0; grid_idx < grids_curr.shape(1); grid_idx++)
            {
                std::array<size_t, 2> px_idx{grids_curr(0, grid_idx), grids_curr(1, grid_idx)};
                mat_a(grid_idx, 0) = sobel_x(px_idx[0], px_idx[1]);
                mat_a(grid_idx, 1) = sobel_y(px_idx[0], px_idx[1]);
                vec_b(grid_idx, 0) = img_curr(px_idx[0], px_idx[1]) - bilinear(grids_next(Col(grid_idx)), img_next)(0,0);
            }

            Matrix<float> solve_a(mat_a.T().matmul(mat_a));
            Matrix<float> solve_b(mat_a.T().matmul(vec_b));

            Matrix<float> increment = qr::solve(solve_a, solve_b);
            if(increment.norm() > prev_inc) { std::cout << "inc inc" << std::endl; break; }
            pt_next += increment;
            // std::cout << "increment: " << increment.norm() << std::endl;
            if(increment.norm() < tol) break;
        }
        pts_next(Col(pt_idx)) = pt_next;
    }
    return pts_next;
}

} // namespace mxm


#endif // _CV_OPTICAL_FLOW_H_
