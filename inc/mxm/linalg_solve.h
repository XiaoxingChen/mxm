#if !defined(_LINALG_SOLVE_H)
#define _LINALG_SOLVE_H

#include "mxm/linalg_mat_ref.h"
#include <iostream>


namespace mxm
{

template<typename DType>
Matrix<DType> solveLUTriangle(const Matrix<DType>& mat, const Matrix<DType>& b, bool l_tri)
{
    if(!mat.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    if(mat.shape(0) != b.shape(0))
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    size_t idx_start = l_tri ? 0 : b.shape(0) - 1;
    int step = l_tri ? 1 : - 1;
    Matrix<DType> x(Matrix<DType>::zeros(b.shape()));

    for(size_t w = 0; w < x.shape(1); w++)
    {
        for(size_t i = idx_start, _1 = 0; _1 < b.shape(0); _1 ++, i += step)
        {
            FloatType off_diag_sum = 0;
            for(size_t j = idx_start, _2 = 0; _2 < _1; _2 ++, j += step)
            {
                off_diag_sum += mat(i, j) * x(j,w);
            }
            x(i,w) = (b(i,w) - off_diag_sum) / mat(i,i);
        }
    }

    return x;
}


template<typename DType>
Matrix<DType> solveLowerTriangle(const Matrix<DType>& lower, const Matrix<DType>& b)
{
    return solveLUTriangle(lower, b, 1);
}

template<typename DType>
Matrix<DType> solveUpperTriangle(const Matrix<DType>& upper, const Matrix<DType>& b)
{
    return solveLUTriangle(upper, b, 0);
}

namespace qr
{
// u: column vector, direction
// a: column vectors
template<typename DType>
Matrix<DType> project(const Matrix<DType>& u, const Matrix<DType>& a)
{
    Matrix<DType> u_dir = u.normalized();
    return u_dir.dot(u_dir.T()).dot(a);
}

template<typename DType>
Matrix<DType> calcMatQFromRotation(const Matrix<DType>& mat_in)
{
    // reference: https://www.math.usm.edu/lambers/mat610/sum10/lecture9.pdf
    if(!mat_in.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    Matrix<DType> mat(mat_in);
    const size_t& n = mat.shape(0);
    Matrix<DType> rot = Matrix<DType>::Identity(n);

    for(size_t j = 0; j < n; j++)
    {
        for(size_t i = n-1; i > j; i--)
        {
            // std::cout << "(i,j): " << i << "," << j << std::endl;
            DType theta = -atan2(mat(i,j), mat(i-1,j));
            // std::cout << "theta: " << theta << std::endl;
            Matrix<DType> sub_rot = Matrix<DType>::Identity(n);

            Matrix<DType> so2({2,2},
                {cos(theta), -sin(theta),
                sin(theta), cos(theta)});

            sub_rot.setBlock(i-1, i-1, so2);
            // std::cout << sub_rot.str() << std::endl;

            rot = sub_rot.matmul(rot);
            mat = sub_rot.matmul(mat);
        }
    }
    return rot.T();
}

template<typename DType>
Matrix<DType> calcMatQFromReflection(const Matrix<DType>& mat)
{
    //reference: https://en.wikipedia.org/wiki/QR_decomposition
    if(!mat.square())
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));

    Matrix<DType> mat_u({mat.shape(0), mat.shape(1)});
    for(size_t i = 0; i < mat.shape(0); i++)
    {
        mat_u(Col(i)) = mat(Col(i));
        for(size_t j = 0; j < i; j++)
        {
            mat_u(Col(i)) -= project(mat_u(Col(j)), mat(Col(i)));
        }
    }

    // std::cout << "Bq: \n" << mat_u.str();
    for(size_t i = 0; i < mat.shape(0); i++)
    {
        mat_u(Col(i)) = mat_u(Col(i)).normalized();
    }

    return mat_u;
}

template<typename DType>
Matrix<DType> calcMatQ(const Matrix<DType>& mat)
{
    return calcMatQFromRotation(mat);
}

template<typename DType>
Matrix<DType> solve(const Matrix<DType>& mat_a, const Matrix<DType>& b)
{
    Matrix<DType> mat_q(calcMatQ(mat_a));

    if((mat_q.matmul(mat_q)).norm() - mat_q.shape(0) > eps())
    {
        // singular
        std::cout << "Singular Matrix!" << std::endl;
        return Matrix<DType>::zeros(b.shape());
    }

    Matrix<DType> mat_r(mat_q.T().matmul(mat_a));
    Matrix<DType> x(solveUpperTriangle(mat_r, mat_q.T().matmul(b)));
    return x;
}
} // namespace qr

template <typename DType>
DType Matrix<DType>::det() const
{
    if(square() && shape(0) == 2)
        return (*this)(0,0) * (*this)(1,1) - (*this)(1,0)*(*this)(0,1);

    Matrix<DType> mat_q(qr::calcMatQ(*this));

    // check full rank
    if((mat_q.matmul(mat_q)).trace() - shape(0) > eps())
        return 0.;

    Matrix<DType> mat_r(mat_q.T().matmul(*this));
    DType det(1);
    for(size_t i = 0; i < shape(0); i++) det *= mat_r(i,i);
    return det;
}

template <typename DType>
Matrix<DType> Matrix<DType>::inv() const
{
    if(!square()) return Matrix<DType>::zeros(shape());
    return qr::solve(*this, Matrix<DType>::Identity(shape(0)));
}

} // namespace mxm



#endif // _LINALG_SOLVE_H
