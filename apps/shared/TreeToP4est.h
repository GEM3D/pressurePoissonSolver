#ifndef TREETOP4EST
#define TREETOP4EST
#include <p4est.h>
#include <OctTree.h>
class TreeToP4est{
    public:
    p4est_connectivity_t* conn;
    p4est_t* p4est;
    TreeToP4est(Tree<2> t);
};
#endif
