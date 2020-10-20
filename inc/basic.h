#ifndef BASIC_H
#define BASIC_H

enum class FLM_BC_MATH : int
{
    Dirichlet = 0,
    Neumann = 1,
    Robin = 2
};

enum class FLM_BC_PHY : int
{
    Wall = 0,
    Inlet = 1,
    Outlet = 2,
    Symmetry = 3
};

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef double FLM_SCALAR;
typedef Eigen::Matrix<FLM_SCALAR, 3, 1> FLM_VECTOR;
typedef Eigen::Matrix<FLM_SCALAR, 3, 3> FLM_TENSOR;

#endif
