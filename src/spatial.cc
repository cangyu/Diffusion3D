#include "../inc/element.h"
#include "../inc/spatial.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

/**
 * Interpolation from cell to node.
 */
void interpolate_nodal_value()
{
    for (auto n : node)
    {
        n->T = 0.0;

        const size_t N = n->cell_dependency.size();
        for (size_t j = 0; j < N; ++j)
        {
            const auto cwf = n->cell_weighting1.at(j);
            auto cdc = n->cell_dependency.at(j);

            n->T += cwf * cdc->T;
        }
    }
}
