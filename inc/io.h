#ifndef IO_H
#define IO_H

#include <istream>
#include <ostream>
#include "basic.h"

void read_mesh(std::istream &fin);

void write_data(std::ostream &out, size_t iter, FLM_SCALAR t);

void read_data(std::istream &in, size_t &iter, FLM_SCALAR &t);

#endif
