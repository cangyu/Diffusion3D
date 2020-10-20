#include "../../inc/element.h"
#include "../../inc/bc.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

void set_bc_desc()
{
    for (auto e : patch)
    {
        if (e->name == "LEFT")
        {
            e->BC = FLM_BC_PHY::Inlet;
            e->T = FLM_BC_MATH::Dirichlet;
        }
        else if (e->name == "RIGHT")
        {
            e->BC = FLM_BC_PHY::Outlet;
            e->T = FLM_BC_MATH::Dirichlet;
        }
        else if (e->name == "WALL")
        {
            e->BC = FLM_BC_PHY::Wall;
            e->T = FLM_BC_MATH::Neumann;
        }
        else
            throw unexpected_patch(e->name);
    }
}

static const FLM_SCALAR T_LEFT = 1500.0; // K
static const FLM_SCALAR T_RIGHT = 300.0; // K

void set_bc_val()
{
    for (auto e : patch)
    {
        if (e->name == "LEFT")
        {
            for (auto f : e->surface)
            {
                f->T = T_LEFT;
            }
        }
        else if (e->name == "RIGHT")
        {
            for (auto f : e->surface)
            {
                f->T = T_RIGHT;
            }
        }
        else if (e->name == "WALL")
        {
            for (auto f : e->surface)
            {
                f->sn_grad_T = 0.0;
            }
        }
        else
            throw unexpected_patch(e->name);
    }
}
