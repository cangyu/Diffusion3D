#include "../inc/element.h"
#include "../inc/green-gauss.h"

extern std::vector<Patch *> patch;
extern std::vector<Node *> node;
extern std::vector<Face *> face;
extern std::vector<Cell *> cell;

/**
 * Calculate gradient on cell centroid.
 * Cell-based.
 * Simple and easy to implement.
 * Only suitable for high-quality mesh.
 */
void gg_m1()
{

}

/**
 * Calculate gradient on cell centroid.
 * Nodal-based.
 * Robust for highly-skewed mesh.
 * Before call to this function, all nodal values should be updated.
 */
void gg_m2()
{
    for(auto f : face)
    {

    }
}