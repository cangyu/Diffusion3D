#include "../inc/element.h"
#include "../inc/error.h"
#include "../inc/lsq.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

typedef Eigen::Matrix<FLM_SCALAR, Eigen::Dynamic, 3> MatX3;
typedef Eigen::Matrix<FLM_SCALAR, 3, Eigen::Dynamic> Mat3X;
typedef Eigen::Matrix<FLM_SCALAR, Eigen::Dynamic, Eigen::Dynamic> MatXX;

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

void calc_lsq_coefficient_matrix()
{
    MatX3 J_rho, J_U[3], J_p, J_p_prime, J_T;

    for (auto c : cell)
    {
        const size_t nF = c->surface.size();

        J_rho.resize(nF, Eigen::NoChange);
        J_U[0].resize(nF, Eigen::NoChange);
        J_U[1].resize(nF, Eigen::NoChange);
        J_U[2].resize(nF, Eigen::NoChange);
        J_p.resize(nF, Eigen::NoChange);
        J_p_prime.resize(nF, Eigen::NoChange);
        J_T.resize(nF, Eigen::NoChange);

        for (size_t j = 0; j < nF; ++j)
        {
            /// Possible coefficients for current face
            auto curFace = c->surface.at(j);
            const auto &d = c->d.at(j);
            const auto w = 1.0 / c->d_mod.at(j);
            if (curFace->at_boundary())
            {
                auto ptc = reinterpret_cast<BoundaryFace *>(curFace)->parent;
                const auto n = c->S.at(j) / curFace->area;

                /// Density
                switch (ptc->rho_BC)
                {
                case FLM_BC_MATH::Dirichlet:
                    J_rho.row(j) << w * d.x(), w * d.y(), w * d.z();
                    break;
                case FLM_BC_MATH::Neumann:
                    J_rho.row(j) << n.x(), n.y(), n.z();
                    break;
                case FLM_BC_MATH::Robin:
                    throw robin_bc_is_not_supported();
                default:
                    break;
                }

                /// Velocity-X
                switch (ptc->U_BC[0])
                {
                case FLM_BC_MATH::Dirichlet:
                    J_U[0].row(j) << w * d.x(), w * d.y(), w * d.z();
                    break;
                case FLM_BC_MATH::Neumann:
                    J_U[0].row(j) << n.x(), n.y(), n.z();
                    break;
                case FLM_BC_MATH::Robin:
                    throw robin_bc_is_not_supported();
                default:
                    break;
                }

                /// Velocity-Y
                switch (ptc->U_BC[1])
                {
                case FLM_BC_MATH::Dirichlet:
                    J_U[1].row(j) << w * d.x(), w * d.y(), w * d.z();
                    break;
                case FLM_BC_MATH::Neumann:
                    J_U[1].row(j) << n.x(), n.y(), n.z();
                    break;
                case FLM_BC_MATH::Robin:
                    throw robin_bc_is_not_supported();
                default:
                    break;
                }

                /// Velocity-Z
                switch (ptc->U_BC[2])
                {
                case FLM_BC_MATH::Dirichlet:
                    J_U[2].row(j) << w * d.x(), w * d.y(), w * d.z();
                    break;
                case FLM_BC_MATH::Neumann:
                    J_U[2].row(j) << n.x(), n.y(), n.z();
                    break;
                case FLM_BC_MATH::Robin:
                    throw robin_bc_is_not_supported();
                default:
                    break;
                }

                /// Pressure
                switch (ptc->p_BC)
                {
                case FLM_BC_MATH::Dirichlet:
                    J_p.row(j) << w * d.x(), w * d.y(), w * d.z();
                    break;
                case FLM_BC_MATH::Neumann:
                    J_p.row(j) << n.x(), n.y(), n.z();
                    break;
                case FLM_BC_MATH::Robin:
                    throw robin_bc_is_not_supported();
                default:
                    break;
                }

                /// Pressure-Correction
                switch (ptc->p_prime_BC)
                {
                case FLM_BC_MATH::Dirichlet:
                    J_p_prime.row(j) << w * d.x(), w * d.y(), w * d.z();
                    break;
                case FLM_BC_MATH::Neumann:
                    J_p_prime.row(j) << n.x(), n.y(), n.z();
                    break;
                case FLM_BC_MATH::Robin:
                    throw robin_bc_is_not_supported();
                default:
                    break;
                }

                /// Temperature
                switch (ptc->T_BC)
                {
                case FLM_BC_MATH::Dirichlet:
                    J_T.row(j) << w * d.x(), w * d.y(), w * d.z();
                    break;
                case FLM_BC_MATH::Neumann:
                    J_T.row(j) << n.x(), n.y(), n.z();
                    break;
                case FLM_BC_MATH::Robin:
                    throw robin_bc_is_not_supported();
                default:
                    break;
                }
            }
            else
            {
                J_rho.row(j) << w * d.x(), w * d.y(), w * d.z(); /// Density
                J_U[0].row(j) << w * d.x(), w * d.y(), w * d.z(); /// Velocity-X
                J_U[1].row(j) << w * d.x(), w * d.y(), w * d.z(); /// Velocity-Y
                J_U[2].row(j) << w * d.x(), w * d.y(), w * d.z(); /// Velocity-Z
                J_p.row(j) << w * d.x(), w * d.y(), w * d.z(); /// Pressure
                J_p_prime.row(j) << w * d.x(), w * d.y(), w * d.z(); /// Pressure-Correction
                J_T.row(j) << w * d.x(), w * d.y(), w * d.z(); /// Temperature
            }
        }

        qr_inv(J_rho, c->J_INV_rho); /// Density
        qr_inv(J_U[0], c->J_INV_U[0]); /// Velocity-X
        qr_inv(J_U[1], c->J_INV_U[1]); /// Velocity-Y
        qr_inv(J_U[2], c->J_INV_U[2]); /// Velocity-Z
        qr_inv(J_p, c->J_INV_p); /// Pressure
        qr_inv(J_p_prime, c->J_INV_p_prime); /// Pressure-Correction
        qr_inv(J_T, c->J_INV_T); /// Temperature
    }
}
