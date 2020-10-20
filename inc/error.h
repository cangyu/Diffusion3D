#ifndef ERROR_H
#define ERROR_H

#include <stdexcept>
#include <string>
#include <cstddef>
#include "basic.h"

struct wrong_face : public std::invalid_argument
{
    explicit wrong_face(size_t idx) :
        std::invalid_argument(std::to_string(idx))
    {}
};

struct unsupported_shape : public std::invalid_argument
{
    unsupported_shape(const std::string &part, size_t idx, int shape) :
        std::invalid_argument("\"" + std::to_string(shape) + "\" on " + part + " " + std::to_string(idx) + ".")
    {}
};

struct invalid_boundary_flag : public std::invalid_argument
{
    invalid_boundary_flag(const std::string &part, size_t idx, int flag) :
        std::invalid_argument("\"" + std::to_string(flag) + "\" on " + part + " " + std::to_string(idx) + ".")
    {}
};

struct failed_to_open_file : public std::runtime_error
{
    explicit failed_to_open_file(const std::string &fn) :
        std::runtime_error(fn)
    {}
};

struct failed_to_create_folder : public std::runtime_error
{
    explicit failed_to_create_folder(const std::string &f) :
        std::runtime_error(f)
    {}
};

struct unsupported_boundary_condition : public std::invalid_argument
{
    explicit unsupported_boundary_condition(FLM_BC_MATH x) :
        std::invalid_argument("math " + std::to_string((int) x))
    {}

    explicit unsupported_boundary_condition(FLM_BC_PHY x) :
        std::invalid_argument("phy " + std::to_string((int) x))
    {}

    explicit unsupported_boundary_condition(const std::string &bc) :
        std::invalid_argument(bc)
    {}
};

struct empty_connectivity : public std::runtime_error
{
    explicit empty_connectivity(int idx) :
        std::runtime_error("Both c0 and c1 are NULL on face " + std::to_string(idx) + ".")
    {}
};

struct inconsistent_connectivity : public std::runtime_error
{
    explicit inconsistent_connectivity(const std::string &msg) :
        std::runtime_error(msg)
    {}
};

struct unexpected_patch : public std::runtime_error
{
    explicit unexpected_patch(const std::string &name) :
        std::runtime_error("\"" + name + "\" is not a pre-defined boundary patch.")
    {}
};

struct no_parent_patch : public std::runtime_error
{
    explicit no_parent_patch(size_t idx) :
        std::runtime_error("NO parent patch for INTERNAL face " + std::to_string(idx) + ".")
    {}
};

struct insufficient_vertexes : public std::runtime_error
{
    explicit insufficient_vertexes(size_t i) :
        std::runtime_error("No enough vertices within cell " + std::to_string(i) + ".")
    {}
};

struct inconsistent_mesh : public std::runtime_error
{
    inconsistent_mesh() :
        std::runtime_error("Input data is NOT consistent with given mesh!")
    {}
};

#endif
