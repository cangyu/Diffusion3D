#ifndef ELEMENT_H
#define ELEMENT_H

#include <cstddef>
#include <vector>
#include <string>
#include "basic.h"
#include "error.h"

class Patch;

class Cell;

class Node
{
public:
    size_t index; /// 1-based global index
    bool at_boundary; /// Boundary flag
    FLM_VECTOR coordinate; /// 3D cartesian location

    /* Connectivity to cells */
    std::vector<Cell *> cell_dependency;
    std::vector<FLM_SCALAR> cell_weighting1; /// 1/||r||
    std::vector<FLM_SCALAR> cell_weighting2; /// 1/||r||^2
    std::vector<FLM_SCALAR> cell_weighting3; /// 1/V

    /* Variable */
    FLM_SCALAR T; /// Temperature
};

class Face
{
public:
    size_t index; /// 1-based global index
    FLM_VECTOR centroid; /// 3D cartesian location of centroid
    FLM_SCALAR area; /// Area of the face

    /* Connectivity to nodes */
    std::vector<Node *> vertex;

    /* Connectivity to cells */
    Cell *cell_dependency[2];
    FLM_SCALAR cell_weighting1[2]; /// 1/||r||
    FLM_SCALAR cell_weighting2[2]; /// 1/||r||^2
    FLM_SCALAR cell_weighting3[2]; /// 1/V
    FLM_VECTOR r[2]; /// Displacement vector
    FLM_VECTOR n[2]; /// Unit normal vector

    /* Property */
    FLM_SCALAR kappa; /// Thermal conductivity

    /* Variable */
    FLM_SCALAR T; /// Temperature

public:
    virtual ~Face() = default;

    virtual bool at_boundary() const = 0;

    virtual void set_parent(Patch *p) = 0;

    Cell *c0() { return cell_dependency[0]; }

    Cell *c1() { return cell_dependency[1]; }

    const FLM_VECTOR &r0() const { return r[0]; }

    const FLM_VECTOR &r1() const { return r[1]; }

    const FLM_VECTOR &n0() const { return n[0]; } /// From c0 to c1

    const FLM_VECTOR &n1() const { return n[1]; } /// From c1 to c0
};

class InternalFace : public Face
{
public:
    /* Gradient */
    FLM_VECTOR grad_T;

public:
    ~InternalFace() = default;

    [[nodiscard]] bool at_boundary() const override { return false; }

    void set_parent(Patch *p) override { throw no_parent_patch(index); }
};

class BoundaryFace : public Face
{
public:
    Patch *parent; /// Connection to high-level

    /* Gradient */
    FLM_SCALAR sn_grad_T; /// In surface outward normal direction

public:
    ~BoundaryFace() = default;

    [[nodiscard]] bool at_boundary() const override { return true; }

    void set_parent(Patch *grp) override { parent = grp; }
};

class Cell
{
public:
    size_t index; /// 1-based global index
    FLM_VECTOR centroid; /// 3D cartesian location of centroid
    FLM_SCALAR volume; /// Volume of the cell

    /* Connectivity to nodes */
    std::vector<Node *> vertex;

    /* Connectivity to faces */
    std::vector<Face *> surface;
    std::vector<FLM_VECTOR> n; /// Surface outward unit normal vector, follow the order in "surface"
    std::vector<FLM_VECTOR> S; /// Surface outward normal vector, follow the order in "surface"
    std::vector<FLM_VECTOR> S_E, S_T; /// Non-Orthogonal decomposition, follow the order in "surface"

    /* Connectivity to cells */
    std::vector<Cell *> cell_adjacency; /// Follow the order in "surface"
    std::vector<FLM_VECTOR> d; /// Displacement vector between adjacent cell centroids
    std::vector<FLM_VECTOR> e; /// Unit displacement vector between adjacent cell centroids
    Eigen::Matrix<FLM_SCALAR, 3, Eigen::Dynamic> J_INV_T; /// Coefficient matrix used by LSQ

    /* Property */
    FLM_SCALAR kappa; /// Thermal conductivity

    /* Variable */
    FLM_SCALAR T; /// Temperature

    /* Gradient */
    FLM_VECTOR grad_T;
};

class Patch
{
public:
    std::string name; /// Identifier

    /* Components */
    std::vector<BoundaryFace *> surface;
    std::vector<Node *> vertex;

    /* B.C. physical classification */
    FLM_BC_PHY BC;

    /* B.C. mathematical specification */
    FLM_BC_MATH T; /// Temperature
};

#endif
