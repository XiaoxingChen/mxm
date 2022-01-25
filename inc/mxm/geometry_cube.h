#if !defined(_GEOMETRY_CUBE_H_)
#define _GEOMETRY_CUBE_H_

#include "linalg.h"

namespace mxm
{
    
template<typename DType=float>
Matrix<DType> generateNCubeVertices(
    const std::vector<DType>& width,
    const std::vector<size_t>& resolutions)
{
    assert(width.size() == resolutions.size());
    size_t dim = resolutions.size();

    std::vector<DType> step(width.size());
    for(size_t i = 0; i < dim; i++) step.at(i) = width.at(i) / resolutions.at(i);

    if(2 == dim)
    {
        std::vector<DType> vertex_data;
        for(size_t i = 0; i < resolutions[0] + 1; i++)
        {
            for(size_t j = 0; j < resolutions[1] + 1; j++)
            {
                vertex_data.push_back(i * step[0]);
                vertex_data.push_back(j * step[1]);     
            }
        }
        return Matrix<DType>(fixRow(2), std::move(vertex_data), COL);
    }
    return Matrix<DType>();
}

template<typename DType=float>
Matrix<size_t> generateNCubeIndices(const std::vector<size_t>& resolutions)
{
    size_t dim = resolutions.size();
    if(2 == dim)
    {
        std::vector<size_t> index_data;
        for(size_t i = 0; i < resolutions[0] + 1; i++)
        {
            for(size_t j = 0; j < resolutions[1] + 1; j++)
            {
                if(j != resolutions[1])
                {
                    index_data.push_back(i * (resolutions[1] + 1) + j);
                    index_data.push_back(i * (resolutions[1] + 1) + j + 1);
                }
                
                if(i != resolutions[0])
                {
                    index_data.push_back(i * (resolutions[1] + 1) + j);
                    index_data.push_back((i + 1) * (resolutions[1] + 1) + j);
                }
                
            }
        }
        return Matrix<size_t>(fixRow(2), std::move(index_data), COL);
    }
    return Matrix<size_t>();
}

} // namespace mxm


#endif // _GEOMETRY_CUBE_H_
