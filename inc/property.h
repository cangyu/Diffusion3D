#ifndef PROPERTY_H
#define PROPERTY_H

#include "basic.h"

FLM_SCALAR Sutherland(FLM_SCALAR T);

void Stokes(FLM_SCALAR mu, const FLM_TENSOR &grad_U, FLM_TENSOR &tau);

#endif
