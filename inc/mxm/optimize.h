#if !defined(__OPTIMIZE_H__)
#define __OPTIMIZE_H__

#include "linalg.h"

namespace mxm
{

namespace opt
{

template<typename DType>
class NoneLinearProblem
{
public:
    DType cost() const { return residual_.T().dot(residual_); }
    const Matrix<DType>& jac() const { return jacobian_; }
    const Vector<DType>& res() const { return residual_; }
    const Vector<DType>& state() const { return state_; }
    void solve(const Vector<DType>& guess, uint8_t step=20, uint8_t verbose=0, std::string method="gn");

    virtual void update(const Vector<DType>& increment) = 0;

protected:
    Vector<DType> state_;
    Matrix<DType> jacobian_;
    Vector<DType> residual_;
};

template<typename DType>
Vector<DType> gaussNewtonIncrement(const NoneLinearProblem<DType>& p, uint8_t verbose)
{
    if(p.jac().norm() < eps())
    {
        std::cout << "WARNING: jac zero!" << std::endl;
        return Vector<DType>::zeros(p.state().size());
    }

    Vector<DType> b = -p.jac().T().matmul(p.res());
    Matrix<DType> hessian = p.jac().T().matmul(p.jac());

    Vector<DType> x = qr::solve(hessian, b);
    return x;
}

template<typename DType>
void NoneLinearProblem<DType>::solve(const Vector<DType>& guess, uint8_t step, uint8_t verbose, std::string method)
{
    state_ = guess;
    Vector<DType> inc;
    for(size_t i = 0; i < step; i++)
    {
        if(method == "gn")
        {
            inc = gaussNewtonIncrement(*this, verbose);
        }
        update(inc);
    }
}

} // namespace opt


} // namespace mxm

#endif // __OPTIMIZE_H__
