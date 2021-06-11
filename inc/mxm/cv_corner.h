#if !defined(_CV_CORNER_H_)
#define _CV_CORNER_H_

#include "cv_pixel.h"
#include "linalg_utils.h"
#include <tuple>

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

        float det = mxm::det(M);
        float tr = M.trace();

        result(i+window_width/2, j+window_width/2) = det - k * tr*tr;
    });
    return result;
}

inline std::array<float, 16> fastCornerCircle(const Matrix<float>& src, size_t idx_x, size_t idx_y)
{
    std::array<float, 16> ret{
        src(idx_x - 3, idx_y + 0),
        src(idx_x - 3, idx_y + 1),
        src(idx_x - 2, idx_y + 2),
        src(idx_x - 1, idx_y + 3),
        src(idx_x - 0, idx_y + 3),
        src(idx_x + 1, idx_y + 3),
        src(idx_x + 2, idx_y + 2),
        src(idx_x + 3, idx_y + 1),
        src(idx_x + 3, idx_y + 0),
        src(idx_x + 3, idx_y - 1),
        src(idx_x + 2, idx_y - 2),
        src(idx_x + 1, idx_y - 3),
        src(idx_x + 0, idx_y - 3),
        src(idx_x - 1, idx_y - 3),
        src(idx_x - 2, idx_y - 2),
        src(idx_x - 3, idx_y - 1)};
    return ret;
}

// Reference:
// https://en.wikipedia.org/wiki/Features_from_accelerated_segment_test
// https://medium.com/data-breach/introduction-to-fast-features-from-accelerated-segment-test-4ed33dde6d65
bool isFastCorner(const std::array<float, 16>& px, float center, float thresh=0.03)
{
    std::vector<size_t> idx_gt;
    std::vector<size_t> idx_lt;
    std::vector<size_t> base_idx{0,4,8,12};
    const size_t CONSECUTIVE_N = 12;
    for(auto i : base_idx)
    {
        if(px.at(i) > center + thresh) idx_gt.push_back(i);
        else if(px.at(i) < center - thresh) idx_lt.push_back(i);
    }
    if(idx_gt.size() < 3 && idx_lt.size() < 3) return false;

    size_t consecutive(0);
    size_t max_consecutive(0);
    int prev=0;
    for(size_t i = 0; i < 2 * px.size(); i++)
    {
        if(px.at(i % px.size()) > center + thresh)
        {
            consecutive = (1 == prev) ? (consecutive + 1) : 0;
            prev = 1;
        }else if(px.at(i % px.size()) < center - thresh)
        {
            consecutive = (-1 == prev) ? (consecutive + 1) : 0;
            prev = -1;
        }else
        {
            consecutive = 0;
            prev = 0;
        }
        if(consecutive >= CONSECUTIVE_N) return true;
    }
    return false;
}

float fastScore(const std::array<float, 16>& pixels, float center)
{
    float score(0);
    for(const auto & px: pixels) score += abs(px - center);
    return score;
}

std::tuple<Matrix<size_t>, Vector<float>>
fastCorners(const Matrix<float>& src, float thresh=0.06)
{
    std::vector<size_t> memory_idx;
    std::vector<float> memory_score;
    for(size_t i = 3; i < src.shape(0) - 3; i++)
    {
        for(size_t j = 3; j < src.shape(1) - 3; j++)
        {
            auto circle_pixels = fastCornerCircle(src, i, j);
            if(!isFastCorner(circle_pixels, src(i,j), thresh)) continue;
            memory_idx.push_back(i);
            memory_idx.push_back(j);
            memory_score.push_back(fastScore(circle_pixels, src(i,j)));
        }
    }

    // return Matrix<size_t>(fixRow(2), std::move(memory_idx), COL);
    return std::make_tuple(
        Matrix<size_t>(fixRow(2), std::move(memory_idx), COL),
        Vector<float>(std::move(memory_score)));
}
} // namespace mxm

#endif // _CV_CORNER_H_
