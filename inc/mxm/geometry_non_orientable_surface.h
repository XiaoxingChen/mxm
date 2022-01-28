#if !defined(_GEOMETRY_NON_ORIENTED_SURFACE_H_)
#define _GEOMETRY_NON_ORIENTED_SURFACE_H_

#include "linalg.h"
#include "rigid_transform.h"

namespace mxm
{

// spin along y axis
// sweep around z axis
template<typename DType=float>
Matrix<DType> generateMobiusStripVertices(
    const DType& width,
    size_t width_resolution,
    const DType& radius,
    size_t radius_resolution)
{
    const size_t dim = 3;
    DType spin_step = M_PI / (radius_resolution);
    DType sweep_step = spin_step * 2;
    DType width_step = width / width_resolution;
    Matrix<DType> point_array({dim, width_resolution + 1});

    for(size_t i = 0; i < width_resolution + 1; i++)
    {
        point_array(0, i) = width_step * i - 0.5 * width;
    }

    Matrix<DType> result({dim, (width_resolution + 1) * radius_resolution});
    RigidTransform<DType, dim> translation({radius,0,0});
    for(size_t i = 0; i < radius_resolution; i++)
    {
        RigidTransform<DType, dim> spin_rot(Rotation<DType, 3>::fromAxisAngle({0,1,0}, spin_step * i));
        RigidTransform<DType, dim> sweep_rot(Rotation<DType, 3>::fromAxisAngle({0,0,1}, sweep_step * i));
        auto local_array = (sweep_rot * translation * spin_rot).apply(point_array);
        result.setBlock(0, (width_resolution + 1) * i, local_array);
    }
    return result;
}

Matrix<size_t> generateMobiusStripLineIndex(
    size_t width_resolution,
    size_t radius_resolution)
{
    std::vector<size_t> index_data;
    index_data.reserve(width_resolution * radius_resolution);
    for(size_t j = 0; j < radius_resolution; j++)
    {
        for(size_t i = 0; i < width_resolution + 1; i++)
        {
            if(i != width_resolution)
            {
                index_data.push_back(j * (width_resolution + 1) + i);
                index_data.push_back(j * (width_resolution + 1) + i + 1);
            }


            if(j+1 == radius_resolution)
            {
                index_data.push_back(j * (width_resolution + 1) + i);
                index_data.push_back(0 * (width_resolution + 1) + (width_resolution - i));
                // index_data.push_back((radius_resolution - j - 1) * (width_resolution + 1) + i);
            }else
            {
                index_data.push_back(j * (width_resolution + 1) + i);
                index_data.push_back((j+1) * (width_resolution + 1) + i);
            }


        }
    }
    return Matrix<size_t>(fixRow(2), index_data, COL);
}

} // namespace mxm


#endif // _GEOMETRY_NON_ORIENTED_SURFACE_H_
