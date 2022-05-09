#if !defined(_OPTIMIZE_SFM_GN_H_)
#define _OPTIMIZE_SFM_GN_H_

#include "linalg.h"

namespace mxm
{
namespace opt
{

template<typename DType>
using BlockDiagMatrix = std::vector<Matrix<DType>>;

template<typename DType>
BlockDiagMatrix<DType>
transpose(const BlockDiagMatrix<DType>& in)
{
    BlockDiagMatrix<DType> ret(in.size());
    for(size_t i = 0; i < in.size(); i++)
    {
        ret.at(i) = in.at(i).T();
    }
    return ret;
}

template<typename DType>
BlockDiagMatrix<DType>
invBlockDiag(const BlockDiagMatrix<DType>& in)
{
    BlockDiagMatrix<DType> ret(in.size());
    for(size_t i = 0; i < in.size(); i++)
    {
        ret.at(i) = inv(in.at(i));
    }
    return ret;
}
template<typename DType>
BlockDiagMatrix<DType>
solveBlockDiag(const BlockDiagMatrix<DType>& mat_a, const Matrix<DType>& b)
{
    BlockDiagMatrix<DType> ret(mat_a.size());
    size_t row_n = 0;
    for(size_t i = 0; i < mat_a.size(); i++)
    {
        ret.at(i) = solve(mat_a.at(i), b({row_n, row_n + mat_a.at(i).shape(0)}, {}));
        row_n += mat_a.at(i).shape(0);
    }
    return ret;
}

template<typename DType> Matrix<DType>
asDense(const BlockDiagMatrix<DType>& block_diag)
{
    Shape ret_shape{0,0};
    for(auto& block:  block_diag)
    {
        ret_shape[0] += block.shape(0);
        ret_shape[1] += block.shape(1);
    }
    Shape up_left{0,0};
    Matrix<DType> ret = Matrix<DType>::zeros(ret_shape);
    for(auto& block:  block_diag)
    {
        ret.setBlock(up_left[0], up_left[1], block);
        up_left[0] += block.shape(0);
        up_left[1] += block.shape(1);
    }
    return ret;
}

// reference: http://ceres-solver.org/nnls_solving.html#dense-schur-sparse-schur
template<typename DType>
Matrix<DType>
solveSfMGaussNewton(
    const Matrix<DType>& jac_dense,
    const BlockDiagMatrix<DType>& jac_sparse,
    const Matrix<DType>& res)
{
    // Vector<DType> b = -p.jac().T().matmul(p.res());
    Matrix<DType> hess_dense = jac_dense.T().matmul(jac_dense);
    Matrix<DType> hess_edge = jac_dense.T().matmul(asDense(jac_sparse));
    BlockDiagMatrix<DType> hess_diag(jac_sparse.size());
    for(size_t i = 0; i < hess_diag.size(); i++)
    {
        hess_diag.at(i) = jac_sparse.at(i).T().matmul(jac_sparse.at(i));
    }
    Matrix<DType> hess_diag_inv = asDense(invBlockDiag(hess_diag));

    Matrix<DType> b_dense = -jac_dense.T().matmul(res);
    Matrix<DType> b_sparse = -asDense(jac_sparse).T().matmul(res);



    Matrix<DType> reduced_cam_mat = hess_dense - hess_edge.matmul( hess_diag_inv ).matmul(hess_edge.T());
    Matrix<DType> reduced_cam_mat_b = b_dense - hess_edge.matmul(hess_diag_inv).matmul(b_sparse);

    Matrix<DType> inc_dense = qr::solve(reduced_cam_mat, reduced_cam_mat_b);
    Matrix<DType> inc_sparse = hess_diag_inv.matmul( b_sparse - hess_edge.T().matmul(inc_dense) );

    Matrix<DType> ret = vstack(inc_dense, inc_sparse);

    if(0){
        Matrix<DType> jac_full = hstack(jac_dense, asDense(jac_sparse) );
        Vector<DType> b = -jac_full.T().matmul(res);
        Matrix<DType> hessian = jac_full.T().matmul(jac_full);

        std::cout << "verify: \n" << mxm::to_string (hessian.matmul(ret)) << std::endl;
        std::cout << "b: \n" << mxm::to_string (b) << std::endl;
    }
    return ret;
}


} // namespace opt

} // namespace mxm


#endif // _OPTIMIZE_SFM_GN_H_
