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
        if (e->name == "UP")
        {
            e->BC = FLM_BC_PHY::Wall;
            e->T = FLM_BC_MATH::Dirichlet;
        }
        else if (e->name == "DOWN")
        {
            e->BC = FLM_BC_PHY::Wall;
            e->T = FLM_BC_MATH::Dirichlet;
        }
        else if (e->name == "LEFT")
        {
            e->BC = FLM_BC_PHY::Wall;
            e->T = FLM_BC_MATH::Neumann;
        }
        else if (e->name == "RIGHT")
        {
            e->BC = FLM_BC_PHY::Wall;
            e->T = FLM_BC_MATH::Neumann;
        }
        else if (e->name == "FRONT")
        {
            e->BC = FLM_BC_PHY::Wall;
            e->T = FLM_BC_MATH::Neumann;
        }
        else if (e->name == "BACK")
        {
            e->BC = FLM_BC_PHY::Wall;
            e->T = FLM_BC_MATH::Neumann;
        }
        else
            throw unexpected_patch(e->name);
    }
}

static const FLM_SCALAR T_UP = 1500.0; // K
static const FLM_SCALAR T_DOWN = 300.0; // K

void set_bc_val()
{
    for (auto e : patch)
    {
        if (e->name == "UP")
        {
            for (auto f : e->surface)
            {
                f->T = T_UP;
            }
        }
        else if (e->name == "DOWN")
        {
            for (auto f : e->surface)
            {
                f->T = T_DOWN;
            }
        }
        else if (e->name == "LEFT")
        {
            for (auto f : e->surface)
            {
                f->sn_grad_T = 0.0;
            }
        }
        else if (e->name == "RIGHT")
        {
            for (auto f : e->surface)
            {
                f->sn_grad_T = 0.0;
            }
        }
        else if (e->name == "FRONT")
        {
            for (auto f : e->surface)
            {
                f->sn_grad_T = 0.0;
            }
        }
        else if (e->name == "BACK")
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
