#if !defined(_INTERPOLATION_H_)
#define _INTERPOLATION_H_

#include "linalg_vec.h"


namespace mxm
{

namespace interp
{

//
// https://en.wikipedia.org/wiki/Bilinear_interpolation#Unit_square
//
// FType: should be either float or double
// DType: can be any type that overloaded with DType::operator*(PType scalar)
//      for example: PixelRGB, Complex, Quaternion
// square: matrix with shape {2,2}
//      | f(0,0), f(0,1) |
//      | f(1,0), f(1,1) |
// pos: (x,y) position
template<typename FType, typename DType>
Matrix<DType>
bilinearUnitSquare(const Matrix<FType>& pos, const Matrix<DType>& square, const std::string& padding="nearest", const DType& padding_val=Traits<DType>::zero() )
{
#if 0
    if(square.shape() != Shape{2,2})
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    Matrix<FType> vx({2,1},{FType(1) - pos(0,0), pos(0,0)});
    Matrix<FType> vy({2,1},{FType(1) - pos(1,0), pos(1,0)});
    return vx.T().matmul(square).matmul(vy)(0,0);
#else
    Matrix<DType> ret({pos.shape(1), 1});
    for(size_t i = 0; i < pos.shape(1); i++)
    {
        float x = pos(0,i);
        float y = pos(1,i);
        if(0 <= x && x + 1 < square.shape(0) && 0 <= y && y + 1 < square.shape(1))
        {
            FType l0(std::floor(x) + 1 - x);
            FType l1(1 - l0);
            FType r0(std::floor(y) + 1 - y);
            FType r1(1 - r0);
            size_t c0(x);
            size_t c1(y);

            ret(i,0) = l0 * square(c0, c1) * r0 + l1 * square(c0+1, c1) * r0
                    + l0 * square(c0, c1+1) * r1 + l1 * square(c0+1, c1+1) * r1;
        }else
        {
            if(std::string("zero") == padding)
            {
                ret(i,0) = Traits<FType>::zero();
            }else if(std::string("nearest") == padding)
            {
                x = std::min<FType>(square.shape(0) - 1, std::max<FType>(0, x));
                y = std::min<FType>(square.shape(1) - 1, std::max<FType>(0, y));
                ret(i,0) = square(size_t(x), size_t(y));
            }else
            {
                ret(i,0) = padding_val;
            }

        }
    }
    return ret;
#endif
}

//
// https://codeplea.com/triangular-interpolation
//
// DType: should be either float or double
// v: vertex matrix with shape {2, 3}, dim is the dimension of space
// val: shape {N, 3}
template<typename DType>
Matrix<DType> triangular(const Vector<DType>& pos, const Matrix<DType>& v, const Matrix<DType>& val)
{
    Vector<DType> weight(3);
    constexpr size_t X = 0;
    constexpr size_t Y = 1;
    DType dom = (v(Y,1) - v(Y,2)) * (v(X,0) - v(X,2)) + (v(X,2) - v(X,1)) * (v(Y,0) - v(Y,2));
    weight(0) = ((v(Y,1) - v(Y,2)) * (pos(X) - v(X,2)) + (v(X,2) - v(X,1)) * (pos(Y) - v(Y,2))) / dom;
    weight(1) = ((v(Y,2) - v(Y,0)) * (pos(X) - v(X,2)) + (v(X,0) - v(X,2)) * (pos(Y) - v(Y,2))) / dom;
    weight(2) = 1 - weight(0) - weight(1);

    return val.matmul(weight);
}


} // namespace interp
} // namespace mxm




#endif // _INTERPOLATION_H_
