#if !defined(_CV_BASIC_H_)
#define _CV_BASIC_H_

#include "cv_pixel.h"
#include "cv_kernel.h"
#include "interpolation.h"
#include "linalg_utils.h"
#include "model_camera.h"

#include <functional>
#include <map>
namespace mxm
{
template<typename PType>
Matrix<PType>
resize(const Matrix<PType>& img, const Shape& shape, const std::string& strategy="bilinear")
{
    Matrix<PType> ret(shape);

    using DType = typename Traits<PType>::ArithType;
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

template<typename PType>
Matrix<PType> flip(const Matrix<PType>& img)
{
    Matrix<PType> ret(img.shape());
    ret.traverse([&](auto i, auto j){
        ret(i,j) = img(i, img.shape(1) - 1 - j);
    });
    return ret;
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
    return Matrix<size_t>(fixRow(2), std::move(local_maximums), COL);
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
    return Matrix<size_t>(fixRow(2), std::move(coords), COL);
}

inline std::vector<std::array<size_t, 2>> greaterThan(
    const Matrix<float>& img, float thresh)
{
    std::vector<std::array<size_t, 2>> coord_buffer;
    img.traverse([&](auto i, auto j){ if(img(i,j) > thresh) coord_buffer.push_back({i,j}); });
    std::sort(coord_buffer.begin(), coord_buffer.end(),
        [&](auto & lhs, auto & rhs) {return img(lhs[0], lhs[1]) < img(rhs[0], rhs[1]);});
    return coord_buffer;
}

inline Matrix<size_t> gridPartitionNonMaximalSuppression(
    const std::vector<std::array<size_t, 2>>& coord_buffer,
    size_t max_pt_num,
    float min_dist)
{
    using Coord2D = std::array<size_t, 2>;

    size_t grid_width = size_t(min_dist + 0.5);
    float dist_square = min_dist * min_dist;
    // grid_coord to idx_of_idx
    std::map<Coord2D, std::vector<Coord2D>> coord_grid_to_img_map;
    auto isOverlap = [&](
        const Coord2D& img_coord,
        const Coord2D& grid_coord,
        float dist_square)->bool
        {
            if(0 == coord_grid_to_img_map.count(grid_coord)) return false;

            for(auto & stored_img_coord : coord_grid_to_img_map[grid_coord])
            {
                float dx = float(img_coord[0]) - float(stored_img_coord[0]);
                float dy = float(img_coord[1]) - float(stored_img_coord[1]);
                if((dx *dx + dy*dy) < dist_square) return true;
            }
            return false;
        };

    std::vector<size_t> result_mem;

    for(auto & img_coord : coord_buffer)
    {
        size_t grid_x = img_coord[0] / grid_width;
        size_t grid_y = img_coord[1] / grid_width;
        Coord2D grid_coord;

        bool has_overlap = false;
        for(grid_coord[0] = grid_x > 0 ? grid_x-1 : 0; grid_coord[0] < grid_x+2; grid_coord[0]++)
        {
            for(grid_coord[1] = grid_y > 0 ? grid_y-1 : 0; grid_coord[1] < grid_y+2; grid_coord[1]++)
            {
                has_overlap = has_overlap || isOverlap(img_coord, grid_coord, dist_square);
                if(has_overlap) break;
            }
            if(has_overlap) break;
        }
        // std::cout << grid_x << "," << grid_y << std::endl;
        if(!has_overlap)
        {
            coord_grid_to_img_map[{grid_x, grid_y}].push_back(img_coord);
            result_mem.push_back(img_coord[0]);
            result_mem.push_back(img_coord[1]);
        }
        else
        {
            // std::cout << "overlap" << std::endl;
        }
        if((result_mem.size() / 2) >= max_pt_num) break;
    }

    return Matrix<size_t>(fixRow(2), std::move(result_mem), COL);
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
    return Matrix<size_t>(fixRow(2), std::move(corners), COL);
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
    return interp::bilinearUnitSquare(pts, img);
#endif
}

template<typename PType> using ImagePyramid = std::vector<Matrix<PType>>;
template<typename PType>
std::vector<Matrix<PType>> gaussPyramid(const Matrix<PType>& img, size_t level, bool return_blured=true)
{
    using DType = typename Traits<PType>::ArithType;
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

inline std::vector<std::ptrdiff_t> halfQuarterBresenhamCircle_(float radius)
{
    std::vector<std::ptrdiff_t> ret;
    std::ptrdiff_t max_y = std::ptrdiff_t(sqrt(0.5f) * radius + 0.5);
    ret.reserve(2 * (max_y * 8));
    float r2 = radius * radius;
    for(std::ptrdiff_t y = 0; y <= max_y; y++)
    {
        std::ptrdiff_t x = sqrt(r2 - y*y) + 0.5;
        if(y == max_y && x == ret.at(ret.size() - 2)) continue;
        ret.push_back(x);
        ret.push_back(y);
    }
    return ret;
}

inline void expendFullBresenhamCircle_(std::vector<ptrdiff_t>& buff)
{
    // to quarter
    int idx = buff.size() - 2;

    while(idx >= 0)
    {
        if(buff.at(idx) != buff.at(idx + 1))
        {
            buff.push_back(buff.at(idx + 1));
            buff.push_back(buff.at(idx + 0));
        }
        idx -= 2;
    }
    // to half
    idx = buff.size() - 2;
    while(idx >= 0)
    {
        if(buff.at(idx) != 0)
        {
            buff.push_back(-buff.at(idx + 0));
            buff.push_back( buff.at(idx + 1));
        }
        idx -= 2;
    }
    // to full
    idx = buff.size() - 2;
    while(idx >= 0)
    {
        if(buff.at(idx + 1) != 0)
        {
            buff.push_back( buff.at(idx + 0));
            buff.push_back(-buff.at(idx + 1));
        }
        idx -= 2;
    }
}

// Reference:
// https://www.geeksforgeeks.org/bresenhams-circle-drawing-algorithm/
inline Matrix<ptrdiff_t> bresenhamCircle(float radius)
{
    auto circle = halfQuarterBresenhamCircle_(radius);
    expendFullBresenhamCircle_(circle);
    return Matrix<ptrdiff_t>(fixRow(2), std::move(circle), COL);
}

// return vector, index as x, value as y.
inline std::vector<size_t> bresenhamHeight_(float radius)
{
    std::vector<size_t> ret;
    size_t max_x = size_t(radius + 0.5);
    float r2 = radius * radius;
    ret.push_back(max_x);


    auto calcErr = [&](float r2_, float x_, float y_) {
        return abs(r2_ - (x_*x_ + y_*y_));
    };
    while(ret.size() < max_x + 1)
    {
        float e0 = calcErr(r2, ret.size(), ret.back());
        float e1 = calcErr(r2, ret.size(), ret.back() - 1);
        ret.push_back(e0 < e1 ? ret.back() : ret.back() - 1);
    }

    return ret;
}

inline void traverseBresenhamCircleArea(
    float radius,
    const std::array<size_t, 2>& center,
    std::function<void(size_t, size_t)> f)
{
    assert(center[0] > radius);
    assert(center[1] > radius);

    auto hs = bresenhamHeight_(radius);
    // std::cout << mxm::to_string(hs) << std::endl;
    for(auto i = 0; i < hs.size(); i++)
    {
        for(int h = -1 * hs.at(i); h < int(hs.at(i) + 1); h++)
        {
            f(center[0] + h, center[1] + i);
            if(0 != i) f(center[0] + h, center[1] - i);
        }
    }
}

template<typename DType>
Matrix<DType> homographicalConstrains(
    const Matrix<DType>& pts_src,
    const Matrix<DType>& pts_dst)
{
    assert(pts_src.shape() == pts_dst.shape());
    const size_t HOMO_SIZE = 9;
    Matrix<DType> ret({2 * pts_src.shape(1), HOMO_SIZE});
    for(size_t i = 0; i < pts_src.shape(1); i++)
    {
        ret(Block({i*2, i*2+2}, {0, end()})) = Matrix<float>({2,HOMO_SIZE},{
            -pts_src(0, i), -pts_src(1, i), -1, 0, 0, 0, pts_src(0, i) * pts_dst(0, i), pts_src(1, i) * pts_dst(0, i), pts_dst(0, i),
            0, 0, 0, -pts_src(0, i), -pts_src(1, i), -1, pts_src(0, i) * pts_dst(1, i), pts_src(1, i) * pts_dst(1, i), pts_dst(1, i)
        }, ROW);
    }

    return ret;
}

// TODO: Smallest eigenvector by power method, see [3].
// Reference:
// [1] https://medium.com/all-things-about-robotics-and-computer-vision/homography-and-how-to-calculate-it-8abf3a13ddc5
// [2] http://www.cs.cmu.edu/~16385/s17/Slides/10.2_2D_Alignment__DLT.pdf
// [3] https://math.stackexchange.com/questions/271864/how-to-compute-the-smallest-eigenvalue-using-the-power-iteration-algorithm

template<typename DType>
Matrix<DType> findHomographyMatrix(
    const Matrix<DType>& pts_src,
    const Matrix<DType>& pts_dst)
{
    auto H = homographicalConstrains(pts_src, pts_dst);
    Matrix<float> v;
    if(pts_src.shape(1) > 4)
    {
        auto eig_val_vec = symmetricEig(H.T().matmul(H));
        v = eig_val_vec[1](Col(end() - 1));
    }else
    {
        v = orthogonalComplement(H).T();
    }

    // std::cout << mxm::to_string(v.T() * (1.f / v(8,0))) << std::endl;
    return Matrix<DType>({3,3},{
        v(0,0),v(1,0),v(2,0),
        v(3,0),v(4,0),v(5,0),
        v(6,0),v(7,0),v(8,0)}) * (DType(1) / v(8,0));
}

#if 0
template<typename PType>
Matrix<PType> undistort(
    const Matrix<PType>& img_src,
    const Camera<typename Traits<PType>::ArithType, 3>& cam,
    const RadialTangentialDistortion<typename Traits<PType>::ArithType>& distortion)
{
    Matrix<PType> img_out(img_src.shape());
    img_out.traverse([&](auto i, auto j){
        auto dir = cam.pixelDirection(Vector<size_t>{i,j});
        auto coord = cam.matrix().matmul(distortion.distort(dir));
        img_out(i,j) = interp::bilinearUnitSquare(coord, img_src, "zero")(0,0);
    });
    return img_out;
}
#endif

} // namespace mxm

#endif // _CV_BASIC_H_
