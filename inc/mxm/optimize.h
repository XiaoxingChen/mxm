#if !defined(__OPTIMIZE_H__)
#define __OPTIMIZE_H__

#include "linalg.h"

namespace mxm
{
namespace opt
{
template<typename DType> class NoneLinearProblem;
} // namespace opt

template<typename DType>
std::string to_string(const opt::NoneLinearProblem<DType>& problem, size_t prec);
namespace opt
{

template<typename DType>
class NoneLinearProblem
{
public:
    DType cost() const { return residual_.T().matmul(residual_)(0,0); }
    const Matrix<DType>& jac() const { return jacobian_; }
    const Vector<DType>& res() const { return residual_; }
    // const Vector<DType>& state() const { return state_; }
    void solve(uint8_t step=20, uint8_t verbose=0, std::string method="gn");

    virtual void update(const Vector<DType>& increment) = 0;

protected:
    // Vector<DType> state_;
    Matrix<DType> jacobian_;
    Vector<DType> residual_;
};

template<typename DType>
Vector<DType> gaussNewtonIncrement(const NoneLinearProblem<DType>& p, uint8_t verbose)
{
    if(p.jac().norm() < eps())
    {
        std::cout << "WARNING: jac zero!" << std::endl;
        return Vector<DType>::zeros(p.res().size());
    }

    Vector<DType> b = -p.jac().T().matmul(p.res());
    Matrix<DType> hessian = p.jac().T().matmul(p.jac());

    Vector<DType> x = qr::solve(hessian, b);
    if(verbose > 1) std::cout << "inc: " << mxm::to_string(x.T());

    if(verbose > 3) std::cout << "b: " << mxm::to_string(b);
    if(verbose > 4) std::cout << "hessian:\n" << mxm::to_string(hessian);
    return x;
}

template<typename DType>
void NoneLinearProblem<DType>::solve(uint8_t step, uint8_t verbose, std::string method)
{
    std::cout << mxm::to_string(*this, verbose);

    Vector<DType> inc;
    for(size_t i = 0; i < step; i++)
    {
        if( cost() < 10 * eps<DType>()) break;
        if(method == "gn")
        {
            inc = gaussNewtonIncrement(*this, verbose);
        }
        update(inc);
        if(verbose > 0) std::cout << "\ni: " << i << std::endl;
        std::cout << mxm::to_string(*this, verbose);
    }
}

} // namespace opt

template<typename DType>
std::string to_string(const opt::NoneLinearProblem<DType>& problem, size_t prec)
{
    std::string ret;
    if(prec > 0) ret += std::string("cost: ") + std::to_string(problem.cost()) + "\n";

    if(prec > 2) ret += std::string("jac:\n") + mxm::to_string(problem.jac());
    if(prec > 3) ret += std::string("res:\n") + mxm::to_string(problem.res());
    return ret;
}

} // namespace mxm

#endif // __OPTIMIZE_H__
