#include "../../inc/element.h"
#include "../../inc/ic.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

static const FLM_SCALAR T0 = 300.0; /// K

void zero_init()
{
    /// Cell
    for (auto c : cell)
    {
        c->T = T0;
    }

    /// Internal Face
    for (auto f : face)
    {
        if (!f->at_boundary())
        {
            f->T = T0;
        }
    }
}
