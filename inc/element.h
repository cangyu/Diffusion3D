#ifndef ELEMENT_H
#define ELEMENT_H

#include <cstddef>
#include <vector>
#include <array>
#include <string>
#include "basic.h"
#include "error.h"

class Patch;

class Cell;

class Node
{
public:
    /// 1-based global index
    size_t index;

    /// Boundary flag
    bool at_boundary;

    /// 3D cartesian location
    FLM_VECTOR coordinate;

    /// Connectivity to cells
    std::vector<Cell *> cell_dependency;

    /// Weighting coefficients for cells
    std::vector<FLM_SCALAR> cell_weighting1; /// 1/||r||
    std::vector<FLM_SCALAR> cell_weighting2; /// 1/||r||^2
    std::vector<FLM_SCALAR> cell_weighting3; /// 1/V
    std::vector<FLM_SCALAR> cell_weighting4; /// Linear preserving

    /// Variable
    FLM_SCALAR T;
};

class Face
{
public:
    /// 1-based global index
    size_t index;

    /// 3D cartesian location of centroid
    FLM_VECTOR centroid;

    /// Area of the face
    FLM_SCALAR area;

    /// Connection to high-level
    Patch *parent;

    /// Connectivity to nodes
    std::vector<Node *> vertex;

    /// Connectivity to cells
    Cell *c0, *c1;

    /// Weighting coefficients for cells
    std::array<FLM_SCALAR, 2> cell_weighting1; /// 1/||r||
    std::array<FLM_SCALAR, 2> cell_weighting2; /// 1/||r||^2
    std::array<FLM_SCALAR, 2> cell_weighting3; /// 1/V

    /// Displacement vector
    FLM_VECTOR r0; /// From centroid of "c0" to face centroid.
    FLM_VECTOR r1; /// From centroid of "c1" to face centroid.

    /// Unit normal vector
    FLM_VECTOR n01; /// From centroid of "c0" to that of "c1".
    FLM_VECTOR n10; /// From centroid of "c1" to that of "c0".

    /// Skewness factor
    FLM_SCALAR alpha;

    /// Property
    FLM_SCALAR kappa; /// Thermal conductivity

    /// Variable
    FLM_SCALAR T;

public:
    virtual ~Face() = default;

    virtual bool at_boundary() const = 0;
};

class InternalFace : public Face
{
public:
    /// Gradient
    FLM_VECTOR grad_T;

public:
    ~InternalFace() = default;

    bool at_boundary() const override { return false; }
};

class BoundaryFace : public Face
{
public:
    /// Gradient
    FLM_SCALAR sn_grad_T; /// In surface OUTWARD normal direction.

public:
    ~BoundaryFace() = default;

    bool at_boundary() const override { return true; }
};

class Cell
{
public:
    /// 1-based global index
    size_t index;

    /// 3D cartesian location of centroid
    FLM_VECTOR centroid;

    /// Volume of the cell
    FLM_SCALAR volume;

    /// Connectivity to nodes
    std::vector<Node *> vertex;

    /// Connectivity to faces
    std::vector<Face *> surface;

    /// Surface OUTWARD normal vector.
    /// Follow the order in "surface".
    std::vector<FLM_VECTOR> S;

    /// Non-Orthogonal decomposition of "S"
    /// "S" = "S_E" + "S_T"
    /// Follow the order in "surface"
    std::vector<FLM_VECTOR> S_E, S_T;

    /// Connectivity to cells
    /// Follow the order in "surface".
    std::vector<Cell *> cell_adjacency;

    /// Displacement vector between adjacent cell centroids
    /// To face centroid if adjacent cell is empty.
    std::vector<FLM_VECTOR> d;

    /// Property
    FLM_SCALAR kappa; /// Thermal conductivity

    /// Variable
    FLM_SCALAR T;

    /// Gradient
    FLM_VECTOR grad_T;
};

class Patch
{
public:
    /// Identifier
    std::string name;

    /// Included faces
    std::vector<BoundaryFace *> surface;

    /// Included nodes
    std::vector<Node *> vertex;

    /// B.C. physical classification
    FLM_BC_PHY BC;

    /// B.C. mathematical specification
    FLM_BC_MATH T;
};

#endif
