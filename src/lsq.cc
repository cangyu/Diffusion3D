#include "../inc/element.h"
#include "../inc/lsq.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

typedef Eigen::Matrix<FLM_SCALAR, Eigen::Dynamic, 3> MatX3;
typedef Eigen::Matrix<FLM_SCALAR, 3, Eigen::Dynamic> Mat3X;
typedef Eigen::Matrix<FLM_SCALAR, Eigen::Dynamic, Eigen::Dynamic> MatXX;

/// Coefficient matrix
static std::vector<Mat3X> J_INV_T;

/**
 * Convert Eigen's intrinsic QR decomposition matrix into R^-1 * Q^T
 * @param J The coefficient matrix to be factorized.
 * @param J_INV The general inverse of input matrix using QR decomposition.
 */
static void qr_inv(const MatX3 &J, Mat3X &J_INV)
{
    auto QR = J.householderQr();
    const MatXX Q = QR.householderQ();
    const MatXX R = QR.matrixQR().triangularView<Eigen::Upper>();

    const MatXX Q0 = Q.block(0, 0, J.rows(), 3);
    const MatXX R0 = R.block<3, 3>(0, 0);

    J_INV = R0.inverse() * Q0.transpose();
}

void prepare_lsq()
{
    MatX3 J_T;

    /// Allocate storage for coefficient matrix
    J_INV_T.resize(cell.size());

    for (auto c : cell)
    {
        const size_t nF = c->surface.size();

        J_T.resize(nF, Eigen::NoChange);

        for (size_t j = 0; j < nF; ++j)
        {
            /// Possible coefficients for current face
            auto curFace = c->surface.at(j);
            const auto &d = c->d.at(j);
            const auto w = 1.0 / c->d.at(j).norm();
            if (curFace->at_boundary())
            {
                auto ptc = curFace->parent;
                const auto n = c->S.at(j) / curFace->area;

                switch (ptc->T) /// Temperature
                {
                case FLM_BC_MATH::Dirichlet:
                    J_T.row(j) << w * d.x(), w * d.y(), w * d.z();
                    break;
                case FLM_BC_MATH::Neumann:
                    J_T.row(j) << n.x(), n.y(), n.z();
                    break;
                default:
                    throw unsupported_boundary_condition(ptc->T);
                }
            }
            else
            {
                J_T.row(j) << w * d.x(), w * d.y(), w * d.z(); /// Temperature
            }
        }

        qr_inv(J_T, J_INV_T.at(c->index-1)); /// Temperature
    }
}

/**
 * Calculate gradient on cell centroid.
 * Before call to this function:
 *   For Dirichlet boundaries:
 *     Values should be updated;
 *     Surface normal gradient are NOT required.
 *   For Neumann boundaries:
 *     Values are NOT required;
 *     Surface normal gradient should be updated.
 */
void lsq()
{
    /// TODO
}
