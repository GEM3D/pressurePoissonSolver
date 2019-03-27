#ifndef TREETOP4EST
#define TREETOP4EST
#include <p4est.h>
#include <p4est_mesh.h>
#include <OctTree.h>
class TreeToP4est{
    public:
    p4est_connectivity_t* conn;
    p4est_t* p4est;
    p4est_ghost_t* ghost;
    p4est_mesh_t* mesh;
    TreeToP4est(Tree<2> t);
};
#endif
