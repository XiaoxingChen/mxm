#if !defined(_CV_3D_H)
#define _CV_3D_H

#include "cv_basic.h"
#include "rigid_transform.h"

namespace mxm
{

template<typename DType>
Rotation<DType>
icpFindRotation(
    const Matrix<DType>& centralized_pts1,
    const Matrix<DType>& centralized_pts2)
{
    auto mat_w = centralized_pts2.matmul(centralized_pts1.T());
    auto u_d_vt = svd(mat_w);
    Matrix<DType> mat_r = u_d_vt[0].matmul(u_d_vt[2]);

    auto rot = Rotation<DType>::fromMatrix(mat_r);
    return rot;
}

#if 1

template<typename DType>
Matrix<DType> invDiagHomogeneous(const Matrix<DType>& homo)
{
    assert(homo.square() || homo.shape(0) + 1 == homo.shape(1));
    Matrix<DType> ret = Matrix<DType>::zeros(homo.shape());
    for(size_t i = 0; i < (homo.square() ? homo.shape(0) - 1 : homo.shape(0)); i++)
    {
        ret(i, i) = DType(1) / homo(i, i);
        ret(i, homo.shape(1) - 1) = -homo(i, homo.shape(1) - 1) / homo(i, i);
    }
    if(homo.square()) ret(ret.shape(0) - 1, ret.shape(1) - 1) = 1;
    return ret;
}

// normalize matrix: T
// p_2^T T_2^T (T_2^{-T} F T_1^{-1}) T_1 p_1 = 0
// A = (T_2^{-T} F T_1^{-1})
// F = T_2^T A T_1
template<typename DeriveType>
Matrix<typename Traits<DeriveType>::EntryType>
findNormalizeMatrix(
    const MatrixBase<DeriveType>& pts,
    typename Traits<DeriveType>::ArithType target_mean_dist= typename Traits<DeriveType>::ArithType(sqrt(2.)))
{
    using DType = typename Traits<DeriveType>::ArithType;
    DType inv_pt_num = DType(1) / pts.shape(1);
    Matrix<DType> weight_center = sum(pts, 1) * inv_pt_num;
    auto residual = pts - weight_center;
    auto dist_sum = sum(residual*residual, 0);
    dist_sum.traverse([&](auto i, auto j) { dist_sum(i,j) = sqrt(dist_sum(i,j)); });

    DType mean_dist = sum(dist_sum) * inv_pt_num;
    DType k = mean_dist / target_mean_dist;
    Matrix<DType> denorm = Matrix<DType>::identity(pts.shape(0) + 1);
    for(size_t i = 0; i < pts.shape(0); i++)
    {
        denorm(i,i) = k;
        denorm(i, pts.shape(0)) = weight_center(i,0);
    }
    return invDiagHomogeneous(denorm);
}

// Solve homogeneous equation: p_2^T F p_1 = 0
//
// Reference:
// [1] https://en.wikipedia.org/wiki/Eight-point_algorithm
// [2] http://www.cse.psu.edu/~rtc12/CSE486/lecture19.pdf
// [3] https://github.com/opencv/opencv/blob/acc576658ad628d46fd4e79c68c1419d438ce716/modules/calib3d/src/fundam.cpp#L672
template<typename DType>
Matrix<DType>
epipolarEightPoint(const Matrix<DType>& pts1, const Matrix<DType>& pts2)
{
    static const size_t N = 8;
    assert(N == pts1.shape(1));
    assert(pts1.shape() == pts2.shape());
#if 1
    Matrix<DType> coeff({N,N+1});
    for(size_t i = 0; i < N; i++)
    {
        coeff(Row(i)) = Matrix<DType>({1, N+1}, {
            pts2(0, i) * pts1(0, i), pts2(0, i) * pts1(1, i), pts2(0, i),
            pts2(1, i) * pts1(0, i), pts2(1, i) * pts1(1, i), pts2(1, i),
            pts1(0, i),              pts1(1, i),              1 });
    }

    // auto solution = orthogonalComplement(coeff);
    // Matrix<DType> solution = svd(coeff)[2](Row(end() - 1));
    Matrix<DType> solution = symmetricEig(coeff.T().matmul(coeff))[1](Col(end() - 1));
    // std::cout << mxm::to_string(coeff.matmul(solution.T()))  << std::endl;
    // std::cout << mxm::to_string(coeff)  << std::endl;
    // std::cout << mxm::to_string(solution)  << std::endl;
    solution.reshape({3,3}, ROW);

    return solution;
#else
    Matrix<DType> sym = Matrix<DType>::zeros({N+1,N+1});
    for(size_t i = 0; i < N; i++)
    {
        auto r = Matrix<DType>({1, N+1}, {
            pts2(0, i) * pts1(0, i), pts2(0, i) * pts1(1, i), pts2(0, i),
            pts2(1, i) * pts1(0, i), pts2(1, i) * pts1(1, i), pts2(1, i),
            pts1(0, i),              pts1(1, i),              1 });

        sym += r.T().matmul(r);
    }
    Matrix<DType> solution = symmetricEig(sym)[1](Col(end() - 1));
    solution.reshape({3,3}, ROW);
    return solution;
#endif
}

template<typename DType>
Matrix<DType> checkEpipolarConstraints(const Matrix<DType>& mat, const Matrix<DType>& pts1, const Matrix<DType>& pts2)
{
    Matrix<DType> ret({pts1.shape(1), 1});
    for(size_t i = 0; i < pts1.shape(1); i++)
    {
        ret(i, 0) = pts2(Col(i)).T().matmul(mat).matmul(pts1(Col(i)))(0,0);
    }
    return ret;
}
#endif



} // namespace mxm


#endif // _CV_3D_H
