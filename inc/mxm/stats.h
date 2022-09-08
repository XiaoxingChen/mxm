#if !defined(__STATS_H__)
#define __STATS_H__

#include "linalg.h"
namespace mxm
{

template<typename DType>
Matrix<DType> findMean(const Matrix<DType>& data)
{
    Matrix<DType> ret = mxm::sum(data, 1);
    ret *= DType(1.)/DType(data.shape(1));
    return ret;
}

template<typename DType>
Matrix<DType> findCovariance(const Matrix<DType>& data, const Matrix<DType>& mean_vec)
{
    Matrix<DType> residual = data - mean_vec;
    return residual.matmul(residual.T()) / DType(data.shape(1) - 1);
}

namespace pdf
{
template <typename DType>
Vector<DType> gaussian(
    const Matrix<DType>& x,
    const Vector<DType>& mean,
    const Matrix<DType>& cov)
{
    assert(cov.square());
    assert(x.shape(0) == mean.size());

    DType cov_det = det(cov);
    if (cov_det < DType(1e-7))
    {
        return Vector<DType>::zeros(mean.size());
    }
    DType cov_inv = inv(cov);

    DType coeff = DType(1.)/(mxm::sqrt(std::pow(2*M_PI, x.shape(0)) * cov_det));
    Matrix<DType> err_vecs = x - mean;
    Vector<DType> ret;
    for(size_t i = 0; i < mean.size(); i++)
    {
        ret(i) = mxm::exp(DType(-0.5) * err_vecs(Col(i)).transpose() * cov_inv * err_vecs(Col(i))) * coeff;
    }
    return ret;
}
} // namespace pdf


} // namespace mxm

#endif // __STATS_H__
