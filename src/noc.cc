#include "../inc/noc.h"

static inline void Minimum(const FLM_VECTOR &d, const FLM_VECTOR &S, FLM_VECTOR &E)
{
    E = (d.dot(S) / (d.dot(d))) * d;
}

static inline void Orthogonal(const FLM_VECTOR &d, const FLM_VECTOR &S, FLM_VECTOR &E)
{
    E = S.norm() / d.norm() * d;
}

static inline void OverRelaxed(const FLM_VECTOR &d, const FLM_VECTOR &S, FLM_VECTOR &E)
{
    E = (S.dot(S) / d.dot(S)) * d;
}

/**
 * Calculate vectors used for NON-ORTHOGONAL correction locally.
 * @param d Local displacement vector.
 * @param S Local surface outward normal vector.
 * @param E Orthogonal part after decomposing "S".
 * @param T Non-Orthogonal part after decomposing "S", satisfying "S = E + T".
 */
void noc_decompose(const FLM_VECTOR &d, const FLM_VECTOR &S, FLM_VECTOR &E, FLM_VECTOR &T)
{
    OverRelaxed(d, S, E);
    T = S - E;
}
