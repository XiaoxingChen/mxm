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

using FastCornerBresehamCircle = std::array<float, 16>;

inline FastCornerBresehamCircle fastCornerCircle(const Matrix<float>& src, size_t idx_x, size_t idx_y)
{
    auto offset = bresenhamCircle(3);
    FastCornerBresehamCircle ret;
    for(size_t i = 0; i < offset.shape(1); i++)
    {
        ret.at(i) = src(idx_x + offset(i, 0), idx_y + offset(i,1));
    }
    return ret;
}

// Reference:
// https://en.wikipedia.org/wiki/Features_from_accelerated_segment_test
// https://medium.com/data-breach/introduction-to-fast-features-from-accelerated-segment-test-4ed33dde6d65
bool isFastCorner(const FastCornerBresehamCircle& px, float center, float thresh=0.03)
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

float fastScore(const FastCornerBresehamCircle& pixels, float center)
{
    float score(0);
    for(const auto & px: pixels) score += abs(px - center);
    return score;
}

std::vector<std::array<size_t, 2>>
fastCorners(const Matrix<float>& src, float thresh=0.06)
{
    using Coord2D = std::array<size_t, 2>;
    std::map<Coord2D, float> scores;
    std::vector<Coord2D> ret;
    for(size_t i = 3; i < src.shape(0) - 3; i++)
    {
        for(size_t j = 3; j < src.shape(1) - 3; j++)
        {
            auto circle_pixels = fastCornerCircle(src, i, j);
            if(!isFastCorner(circle_pixels, src(i,j), thresh)) continue;
            scores[{i,j}] = fastScore(circle_pixels, src(i,j));
            ret.push_back({i,j});
        }
    }

    std::sort(ret.begin(), ret.end(), [&](auto & lhs, auto& rhs){
        return scores[lhs] > scores[rhs];
    });

    return ret;
}

template<uint8_t K>
const Matrix<int>& briefPatternII();

template<>
const Matrix<int>& briefPatternII<8>()
{
    static const Matrix<int> pattern({4,64},
    {-11,   0, -10,  10,
    -1,   4,  -7,  -2,
    0, -11,  -3,  -7,
    4,   3,  -3,   3,
    5,  -5,  -6,   1,
    1,  -3,  -4,   5,
    2,   1,   3,  -7,
    -3,   8,  -7,  11,
    -2,  -3,   2,   2,
    0,   1,   6,  -4,
    -5,   5, -10,   3,
    4,  -2,   4,  -1,
    10,  -2,   0,   3,
    1,  -3,  -8,  -6,
    6,   8,  -3,  -1,
    -3, -10,  -1,  -1,
    0,  10,  -2,   0,
    -1,  -8,  -6,   5,
    6,   3,   7,   8,
    -6,   4,   5,  -4,
    9,   5,  -4,  10,
    6,   3,  -1,   3,
    -8,  -5,   2,   1,
    -2,  -7,   4,   7,
    5,  -7,   2,  -1,
    -5,  -9,  -7,   9,
    5,   8,  -2,  12,
    1,   3,   7,  -1,
    14,  -1,   0,  -1,
    -2,  15,  -1,  -3,
    -4,  -1,   0,   5,
    15,   5,  -3,   5,
    3,  -2,   0,   1,
    8,  10,   0,   8,
    2,   0,   0,  -3,
    -2,   7,  -9,   0,
    -9,   2,   9,  -7,
    0,   1,   3,  -1,
    -2,  -1,  -4, -12,
    6,   3,  -1,   7,
    -1,  -8,   8,   2,
    10,   0,   3,  -1,
    2,  -4,   3,   0,
    12,   6,  -5,  -6,
    -5,   1,  -4,  -1,
    -1,  16,  -3,  10,
    4,   4,   1,   9,
    0,   6,  -2,   2,
    2,  -9,  -6,   3,
    -9,   1,   1,  -2,
    6,   0,  12,   3,
    4,   1,   4,   2,
    -2,   6,   0,   2,
    -7,  -1,  -9,   2,
    6,   0,  -3,  -4,
    0,   0,  -5,  -8,
    -4,  -5,   3,  12,
    3,  -4,  -3,  -1,
    -1,  10,  -4,   0,
    -9,   0,   5,   0,
    -2,  -9,   1,  11,
    10,  -5,   9,   1,
    1,   5,   2,   0,
    0,  -5,  -4, -13}, COL);
   return pattern;
}

// Reference:
// https://papers.nips.cc/paper/2012/file/59b90e1005a220e2ebc542eb9d950b1e-Paper.pdf
template<uint8_t K>
class BriefDescriptor
{
public:
    using ThisType = BriefDescriptor<K>;
    const uint8_t BYTES = K;
    // BriefDescriptor() { for(size_t i = 0; i < K; i++) data_.at(i) = 0; }
    void operator = (const ThisType& rhs) { data_ = rhs.data_; }

    void setBit(size_t b, bool val)
    {
        uint8_t idx = b / 8;
        uint8_t pos = b % 8;
        data_.at(idx) &= ~(uint8_t(1) << pos);
        data_.at(idx) |= (val << pos);
    }
    uint8_t at(size_t i) const {return data_.at(i);}

    size_t distance(const ThisType& rhs)
    {
        size_t cnt = 0;
        for(size_t i = 0; i < K; i++)
        {
            uint8_t v = at(i) ^ rhs.at(i);
            while(v)
            {
                cnt += (v & 1);
                v >>= 1;
            }
        }
        return cnt;
    }
private:
    std::array<uint8_t, K> data_;
};

template<uint8_t K>
std::string to_string(const BriefDescriptor<K>& v)
{
    std::stringstream stream;
    stream << std::setfill ('0') << std::uppercase;
    stream << std::hex;
    for(size_t i = 0; i < K; i++)
    {
        stream << std::setw(2) << int(v.at(K - 1 - i));
    }
    return stream.str();
}


template<uint8_t K=8>
Matrix<BriefDescriptor<K>> calculateBriefDescriptor(
    const Matrix<float>& img, const Matrix<size_t>& pts)
{
    const size_t patch_size = 32;
    const size_t half_w = patch_size / 2;

    Matrix<BriefDescriptor<K>> ret({pts.shape(1), 1});
    const auto & pattern = briefPatternII<K>();
    for(size_t pt_idx = 0; pt_idx < pts.shape(1); pt_idx++)
    {
        if(pts(0, pt_idx) < half_w
        || pts(0, pt_idx) + half_w >= img.shape(0)
        || pts(1, pt_idx) < half_w
        || pts(1, pt_idx) + half_w >= img.shape(1))
        {
            continue;
        }

        for(size_t b = 0; b < 8 * K; b++)
        {
            bool sign =
                img(pts(0, pt_idx) + pattern(b, 0), pts(1, pt_idx) + pattern(b, 1))
                > img(pts(0, pt_idx) + pattern(b, 2), pts(1, pt_idx) + pattern(b, 3));
            ret(pt_idx, 0).setBit(b, sign);
        }
    }
    return ret;
}


} // namespace mxm

#endif // _CV_CORNER_H_
