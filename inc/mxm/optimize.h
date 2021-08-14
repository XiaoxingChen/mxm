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
    virtual Vector<DType> residualOnPerturb(const Vector<DType>& increment) const { assert(false); return residual_; }
    DType reductionRate(const Vector<DType>& increment) const
    {
        DType res_norm = residual_.norm() ;
        DType actual = residualOnPerturb(increment).norm();
        DType predict = (residual_ + jac().matmul(increment)).norm();
        return (res_norm - actual * actual) / (res_norm - predict * predict);
    }

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

// Reference:
// [1] The Levenberg-Marquardt Algorithm: Implementation and Theory.
//      Link: https://scholar.google.com/scholar_url?url=https://www.osti.gov/servlets/purl/7256021/&hl=en&sa=X&ei=TcsUYeflH8iWywSDwbpQ&scisig=AAGBfm2OFFcmhjfs1cUBu2CvjvQIU6bNRg&oi=scholarr
// [2] https://github.com/RobotLocomotion/eigen-mirror/tree/master/unsupported/Eigen/src/LevenbergMarquardt
template<typename DType>
Vector<DType> levenbergMarquardtIncrement(const NoneLinearProblem<DType>& p, uint8_t verbose)
{
    if(p.jac().norm() < eps<DType>())
    {
        std::cout << "WARNING: jac zero!" << std::endl;
        return Vector<DType>::zeros(p.res().size());
    }

    DType trust_radius = 2.;
    DType damp = 0.05;
    size_t max_it = 1;
    Vector<DType> x = Vector<DType>::zeros(p.jac().shape(1));
    Vector<DType> b = -p.jac().T().matmul(p.res());
    Matrix<DType> hessian = p.jac().T().matmul(p.jac());
    for( size_t i = 0; i < max_it; i++)
    {
        Matrix<DType> hessian_pos = hessian + damp * Matrix<DType>::identity(hessian.shape(0));

        x = qr::solve(hessian_pos, b);


    }

    if(verbose > 1) std::cout << "inc: " << mxm::to_string(x.T());
    if(verbose > 1) std::cout << "damp: " << damp << std::endl;

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
        }else if(method == "lm")
        {
            inc = levenbergMarquardtIncrement(*this, verbose);
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
