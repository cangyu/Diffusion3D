#include "../inc/property.h"

/* Flow condition */
FLM_SCALAR Re = 100.0;

/**
 * Dynamic viscosity of ideal gas.
 * @param T Temperature in Kelvin.
 * @return Dynamic viscosity with unit "Kg / (m * s)".
 */
FLM_SCALAR Sutherland(FLM_SCALAR T)
{
    return 1.45e-6 * std::pow(T, 1.5) / (T + 110.0);
}

/**
 * Viscous shear stress of newtonian fluid under stokes hypothesis.
 * @param mu Dynamic viscosity.
 * @param grad_U Velocity gradient.
 * @param tau Shear stress.
 */
void Stokes(FLM_SCALAR mu, const FLM_TENSOR &grad_U, FLM_TENSOR &tau)
{
    const FLM_SCALAR loc_div3 = (grad_U(0, 0) + grad_U(1, 1) + grad_U(2, 2)) / 3.0;

    tau(0, 0) = 2.0 * mu * (grad_U(0, 0) - loc_div3);
    tau(1, 1) = 2.0 * mu * (grad_U(1, 1) - loc_div3);
    tau(2, 2) = 2.0 * mu * (grad_U(2, 2) - loc_div3);

    tau(0, 1) = tau(1, 0) = mu * (grad_U(0, 1) + grad_U(1, 0));
    tau(1, 2) = tau(2, 1) = mu * (grad_U(1, 2) + grad_U(2, 1));
    tau(2, 0) = tau(0, 2) = mu * (grad_U(2, 0) + grad_U(0, 2));
}
