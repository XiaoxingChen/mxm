#if !defined(__GEOMETRY_TORUS_H__)
#define __GEOMETRY_TORUS_H__

#include "linalg.h"
#include "rotation.h"

namespace mxm
{
template<size_t DIM>
Matrix<size_t> triangulateQuad(const Matrix<size_t>& indices)
{
    assert(false);
    return indices;
}

template<>
Matrix<size_t> triangulateQuad<1>(const Matrix<size_t>& indices)
{
    return indices;
}

template<>
Matrix<size_t> triangulateQuad<2>(const Matrix<size_t>& indices)
{
    assert(indices.shape(0) == 4);
    const size_t TRIANGLE_PER_QUAD = 2;
    Matrix<size_t> ret({3, indices.shape(1) * TRIANGLE_PER_QUAD});
    for(size_t i = 0; i < indices.shape(1); i++)
    {
        ret(0, TRIANGLE_PER_QUAD*i) = indices(0, i);
        ret(1, TRIANGLE_PER_QUAD*i) = indices(1, i);
        ret(2, TRIANGLE_PER_QUAD*i) = indices(2, i);

        ret(0, TRIANGLE_PER_QUAD*i+1) = indices(0, i);
        ret(1, TRIANGLE_PER_QUAD*i+1) = indices(3, i);
        ret(2, TRIANGLE_PER_QUAD*i+1) = indices(2, i);
    }
    return ret;
}

template<>
Matrix<size_t> triangulateQuad<3>(const Matrix<size_t>& indices)
{
    assert(indices.shape(0) == 8);
    const size_t TETRAHEDRON_PER_CUBE = 6;
    Matrix<size_t> ret({4, indices.shape(1) * TETRAHEDRON_PER_CUBE});
    for(size_t i = 0; i < indices.shape(1); i++)
    {
        Vector<size_t> v = indices(Col(i));
        ret.setBlock(0, i*TETRAHEDRON_PER_CUBE,
            Matrix<size_t>({4, TETRAHEDRON_PER_CUBE}, {
                0,1,3,4,
                1,4,5,7,
                1,3,4,7,
                0,2,3,4,
                2,4,6,7,
                0,2,3,7}, COL)
        );
    }
    return ret;
}

template<size_t DIM, typename DType=float>
Matrix<DType> generateNTorus(
    const std::vector<DType>& radius,
    const std::vector<size_t>& resolutions)
{
    assert(radius.size() + 1 == DIM);
    assert(resolutions.size() + 1 == DIM);

    size_t dim = radius.size() + 1;
    std::vector<size_t> vertex_num_by_layer{1};
    for(const auto & r: resolutions)
    {
        vertex_num_by_layer.push_back(vertex_num_by_layer.back() * r);
    }

    Matrix<DType> vertex_buffer = Matrix<DType>({dim, vertex_num_by_layer.back()});
    for(size_t axis = 0; axis < dim - 1; axis++)
    {
        DType angle_per_step = 2 * M_PI / resolutions.at(axis);
        auto plane_axis_0 = Vector<DType>::zeros(dim);
        auto plane_axis_1 = Vector<DType>::zeros(dim);
        Vector<DType> radius_vec(dim);
        plane_axis_0(axis) = 1;
        plane_axis_1(axis+1) = 1;
        radius_vec = plane_axis_0 * radius.at(axis);

        auto src_vertices = vertex_buffer(Block({}, {0, vertex_num_by_layer.at(axis)}));
        src_vertices += radius_vec;
        for(size_t step = 1; step < resolutions.at(axis); step++)
        {
            auto rot = Rotation<DType, DIM>::fromPlaneAngle(plane_axis_0, plane_axis_1, angle_per_step * step);
            auto dst_vertices = vertex_buffer(Block({}, {step * vertex_num_by_layer.at(axis), (step + 1) * vertex_num_by_layer.at(axis)}));
            dst_vertices = rot.apply(src_vertices);
        }
    }
    return vertex_buffer;
}

Matrix<size_t> generateNTorusGraph(const std::vector<size_t>& resolutions)
{
    size_t dim = resolutions.size() + 1;
    std::vector<size_t> vertex_num_by_layer{1};
    for(const auto & r: resolutions)
    {
        vertex_num_by_layer.push_back(vertex_num_by_layer.back() * r);
    }
    Matrix<size_t> quad_index_buffer({1ul << dim, vertex_num_by_layer.back()});
    if(dim == 2) return triangulateQuad<2>(quad_index_buffer);
    else if(dim == 3) return triangulateQuad<3>(quad_index_buffer);
    else if(dim == 4) return triangulateQuad<4>(quad_index_buffer);
    return Matrix<size_t>();
}

Matrix<size_t> generateNTorusToNCubeIndex(const std::vector<size_t>& resolutions)
{
    size_t dim = resolutions.size() + 1;
    std::vector<size_t> index_data;
    if(2 == dim)
    {
        size_t VERTEX_PER_LINE = 2;
        index_data.reserve(resolutions.at(0) * VERTEX_PER_LINE);
        for(size_t i = 0; i < resolutions.at(0); i++)
        {
            size_t i_1 = (i+1 == resolutions.at(0)) ? 0 : i+1;
            index_data.push_back(i);
            index_data.push_back(i_1);
        }
        return Matrix<size_t>({VERTEX_PER_LINE, resolutions.at(0)}, std::move(index_data), COL);
    }else if(3 == dim)
    {
        size_t VERTEX_PER_QUAD = 4;
        size_t quad_num = resolutions.at(0) * resolutions.at(1) * 2;
        index_data.reserve(quad_num * VERTEX_PER_QUAD);
        for(size_t j = 0; j < resolutions.at(1); j++)
        {
            for(size_t i = 0; i < resolutions.at(0); i++)
            {
                size_t i_1 = (i+1 == resolutions.at(0)) ? 0 : i+1;
                size_t j_1 = (j+1 == resolutions.at(1)) ? 0 : j+1;
                index_data.push_back(i + j * resolutions.at(0));
                index_data.push_back(i_1 + j * resolutions.at(0));
                index_data.push_back(i_1 + j_1 * resolutions.at(0));
                index_data.push_back(i + j_1 * resolutions.at(0));
            }
        }
        return Matrix<size_t>({VERTEX_PER_QUAD, quad_num}, std::move(index_data), COL);
    }else if(4 == dim)
    {
        size_t VERTEX_PER_CUBE = 8;
        size_t cube_num = resolutions.at(0) * resolutions.at(1) * resolutions.at(2) * 2;
        index_data.reserve(cube_num * VERTEX_PER_CUBE);
        size_t j_stride = resolutions.at(0);
        size_t k_stride = resolutions.at(0) * resolutions.at(1);
        for(size_t k = 0; k < resolutions.at(2); k++)
        {
            for(size_t j = 0; j < resolutions.at(1); j++)
            {
                for(size_t i = 0; i < resolutions.at(0); i++)
                {
                    size_t i_1 = (i+1 == resolutions.at(0)) ? 0 : i+1;
                    size_t j_1 = (j+1 == resolutions.at(1)) ? 0 : j+1;
                    size_t k_1 = (k+1 == resolutions.at(2)) ? 0 : k+1;

                    index_data.push_back(i + j * j_stride + k * k_stride);
                    index_data.push_back(i_1 + j * j_stride + k * k_stride);
                    index_data.push_back(i_1 + j_1 * j_stride + k * k_stride);
                    index_data.push_back(i + j_1 * j_stride + k * k_stride);
                    index_data.push_back(i + j * j_stride + k_1 * k_stride);
                    index_data.push_back(i_1 + j * j_stride + k_1 * k_stride);
                    index_data.push_back(i_1 + j_1 * j_stride + k_1 * k_stride);
                    index_data.push_back(i + j_1 * j_stride + k_1 * k_stride);
                }
            }
        }
        return Matrix<size_t>({VERTEX_PER_CUBE, cube_num}, std::move(index_data), COL);
    }
    return Matrix<size_t>();
}

Matrix<size_t> generateNTorusLineIndex(const std::vector<size_t>& resolutions)
{
    size_t dim = resolutions.size() + 1;
    std::vector<size_t> index_data;
    const size_t VERTEX_PER_LINE = 2;
    if(2 == dim)
    {
        index_data.reserve(resolutions.at(0) * VERTEX_PER_LINE);
        for(size_t i = 0; i < resolutions.at(0); i++)
        {
            size_t i_1 = (i+1 == resolutions.at(0)) ? 0 : i+1;
            index_data.push_back(i);
            index_data.push_back(i_1);
        }
        return Matrix<size_t>({VERTEX_PER_LINE, resolutions.at(0)}, std::move(index_data), COL);
    }else if(3 == dim)
    {
        size_t line_num = resolutions.at(0) * resolutions.at(1) * 2;
        index_data.reserve(line_num * VERTEX_PER_LINE);
        for(size_t j = 0; j < resolutions.at(1); j++)
        {
            for(size_t i = 0; i < resolutions.at(0); i++)
            {
                size_t i_1 = (i+1 == resolutions.at(0)) ? 0 : i+1;
                size_t j_1 = (j+1 == resolutions.at(1)) ? 0 : j+1;
                index_data.push_back(i + j * resolutions.at(0));
                index_data.push_back(i_1 + j * resolutions.at(0));

                index_data.push_back(i + j * resolutions.at(0));
                index_data.push_back(i + j_1 * resolutions.at(0));
            }
        }
        return Matrix<size_t>({VERTEX_PER_LINE, line_num}, std::move(index_data), COL);
    }else if(4 == dim)
    {
        size_t line_num = resolutions.at(0) * resolutions.at(1) * resolutions.at(2) * 3;
        index_data.reserve(line_num * VERTEX_PER_LINE);
        size_t j_stride = resolutions.at(0);
        size_t k_stride = resolutions.at(0) * resolutions.at(1);
        for(size_t k = 0; k < resolutions.at(2); k++)
        {
            for(size_t j = 0; j < resolutions.at(1); j++)
            {
                for(size_t i = 0; i < resolutions.at(0); i++)
                {
                    size_t i_1 = (i+1 == resolutions.at(0)) ? 0 : i+1;
                    size_t j_1 = (j+1 == resolutions.at(1)) ? 0 : j+1;
                    size_t k_1 = (k+1 == resolutions.at(2)) ? 0 : k+1;

                    index_data.push_back(i + j * j_stride + k * k_stride);
                    index_data.push_back(i_1 + j * j_stride + k * k_stride);

                    index_data.push_back(i + j * j_stride + k * k_stride);
                    index_data.push_back(i + j_1 * j_stride + k * k_stride);

                    index_data.push_back(i + j * j_stride + k * k_stride);
                    index_data.push_back(i + j * j_stride + k_1 * k_stride);
                }
            }
        }
        return Matrix<size_t>({VERTEX_PER_LINE, line_num}, std::move(index_data), COL);
    }
    return Matrix<size_t>();
}
} // namespace mxm

#endif // __GEOMETRY_TORUS_H__
