#if !defined(_ESTIMATOR_KALMAN_FILTER_)
#define _ESTIMATOR_KALMAN_FILTER_

#include "linalg.h"

namespace mxm
{

namespace kf
{

template<typename DType>
std::array<Matrix<DType>, 2>
predict(
    const Matrix<DType>& mean_prev,
    const Matrix<DType>& cov_prev,
    const Matrix<DType>& update,
    const Matrix<DType>& noise,
    const Matrix<DType>& control)
{
    const size_t STATE_DIM(mean_prev.shape(0));
    assert(cov_prev.square()));
    assert(STATE_DIM == cov_prev.shape(0));

    std::array<Matrix<DType>, 2> ret;

    Matrix<DType>& mean_curr(ret[0]);
    Matrix<DType>& cov_curr(ret[1];

    mean_curr = update.matmul(mean_prev) + control;
    cov_curr = update.matmul(cov_prev).matmul(update.T()) + noise;

    return ret;
}

template<typename DType>
Matrix<DType> gain(
    const Matrix<DType>& cov_curr,
    const Matrix<DType>& measure_mat,
    const Matrix<DType>& cov_measure)
{
    Matrix<DType> term2 = mxm::inv(measure_mat.matmul(cov_curr).matmul(measure_mat.T()) + cov_measure);
    Matrix<DType> k = cov_curr.matmul(measure_mat.T()).matmul(term2);
    return k;
}

} // namespace kf


} // namespace mxm


#endif // _ESTIMATOR_KALMAN_FILTER_
