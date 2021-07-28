#if !defined(_CV_3D_H)
#define _CV_3D_H

#include "cv_basic.h"
#include "rigid_transform.h"
#include "optimize.h"
#include "lie_special_orthogonal.h"

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

template<typename DType>
std::array<Matrix<DType>, 2>
basisForUpdate(const Matrix<DType>& v)
{
    Vector<DType> ref_v{1,0,0};
    if(DType(1) - norm(v.T().matmul(ref_v)) < DType(0.01)) ref_v = Vector<DType>{0,1,0};

    auto wedge_v = so::wedge(v);
    Vector<DType> vx = (wedge_v.matmul(ref_v)).normalized();
    Vector<DType> vy = (wedge_v.matmul(vx)).normalized();

    std::array<Matrix<DType>, 2> ret{vx, vy};
    return ret;
}

// The epipolar constraits was modified from p2^T (t^wedge R) p1
// to -p2^T (Rp1)^wedge t
template<typename DType>
class EpipolarOptimize: public opt::NoneLinearProblem<DType>
{
public:
    EpipolarOptimize(const Matrix<DType>& pts1, const Matrix<DType>& pts2)
    :pts1_(pts1), pts2_(pts2), pt_num_(pts1.shape(1))
    {
        assert(pts1.shape() == pts2.shape());
        assert(pts1.shape(1) >= 5);

        this->residual_ = Vector<DType>(pt_num_);
        this->jacobian_ = Matrix<DType>({pt_num_, 5});
    }
    static std::array<Matrix<DType>, 2> rotationFromIncrement(const Matrix<DType>& increment, const Matrix<DType>& t_state);
    static Vector<DType> vecFromTR(const Matrix<DType>& t, const Matrix<DType>& rot);
    // state_: tx ty rx ry rz
    // update: rotate vector t, rotate matrix r
    virtual void update(const Vector<DType>& increment) override
    {
        auto t_r_inc = rotationFromIncrement(increment, t_state_);
        t_state_ = t_r_inc[0].matmul(t_state_);
        r_state_ = t_r_inc[1].matmul(r_state_);

        DType t_norm = mxm::norm(t_state_);
        if(t_norm < eps<DType>())
            t_state_ = Vector<DType>{1,0,0};
        else
            t_state_ /= t_norm;
        t_state_ *= DType(10);

        auto t_basis = basisForUpdate(t_state_);

        for(size_t i = 0; i < pt_num_; i++)
        {
            auto rot_p1_wedge = so::wedge(r_state_.matmul(pts1_(Col(i))));
            auto residual_term = pts2_(Col(i)).T().matmul(rot_p1_wedge);
            // update res
            this->residual_(i) = -residual_term.matmul(t_state_)(0,0);
            // update jac t
            this->jacobian_(i, 0) = residual_term.matmul(so::wedge(t_state_)).matmul((t_basis[0]))(0,0);
            this->jacobian_(i, 1) = residual_term.matmul(so::wedge(t_state_)).matmul((t_basis[1]))(0,0);
            // update jac r
            this->jacobian_(Block({i, i+1}, {2, 5})) = -pts2_(Col(i)).T().matmul(so::wedge(t_state_)).matmul(rot_p1_wedge);
        }
    }

    void initialGuess(const Matrix<DType>& t_init, const Matrix<DType>& r_init)
    {
        t_state_ = t_init;
        r_state_ = r_init;
    }
    const Vector<DType>& tState() const { return t_state_; }
    const Matrix<DType>& rState() const { return r_state_; }

private:
    Matrix<DType> pts1_;
    Matrix<DType> pts2_;
    size_t pt_num_;
    Vector<DType> t_state_;
    Matrix<DType> r_state_;
};

template<typename DType> std::array<Matrix<DType>, 2>
EpipolarOptimize<DType>::rotationFromIncrement(const Matrix<DType>& vec, const Matrix<DType>& t_state)
{
    std::array<Matrix<DType>, 2> t_r;
    auto t_basis = basisForUpdate(t_state);
    t_r[0] = so::exp<3>(
        so::wedge(t_basis[0] * vec(0, 0) + t_basis[1] * vec(1, 0)));
    t_r[1] = so::exp<3>(so::wedge(vec(Block({2,5},{}))));
    return t_r;
}
#if 0
template<typename DType> Vector<DType>
EpipolarOptimize<DType>::vecFromTR(const Matrix<DType>& t, const Matrix<DType>& rot)
{
    Vector<DType> ret(6);
    auto retvec = so::vee(SO::log<3>(rot));
    for(size_t i = 0;i < t.shape(0); i++)
    {
        ret(i) = t(i, 0);
        ret(i + 3) = retvec(i);
    }
    return ret;
}
#endif

template<typename DType>
std::array<Matrix<DType>, 2>
epipolarLeastSquare(const Matrix<DType>& pts1, const Matrix<DType>& pts2, const Vector<DType>& t_guess)
{
    EpipolarOptimize<DType> problem(pts1, pts2);
    problem.initialGuess(t_guess, Matrix<DType>::identity(3));
    problem.solve(10, 0);
    std::array<Matrix<DType>, 2> ret{problem.tState().normalized(), problem.rState()};
    return ret;
    // return EpipolarOptimize<DType>::rotationFromIncrement(problem.state());
}



} // namespace mxm


#endif // _CV_3D_H
