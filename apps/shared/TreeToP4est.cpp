#include "TreeToP4est.h"
#include <p4est_connectivity.h>
#include <p4est_extended.h>
int refine_fn(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quadrant)
{
	Tree<2> &t = *(Tree<2> *) p4est->user_pointer;
	// double x = p4est_root_length
	Node<2> node = t.nodes[t.root];
	for (int i = 1; i <= quadrant->level; i++) {
		if (node.hasChildren()) {
			int quad = (quadrant->x / P4EST_QUADRANT_LEN(i)) & 0x1;
			quad |= ((quadrant->y / P4EST_QUADRANT_LEN(i)) & 0x1) << 1;
			node = t.nodes[node.child_id[quad]];
		} else {
			break;
		}
	}
	return node.hasChildren();
}
TreeToP4est::TreeToP4est(Tree<2> t)
{
	conn  = p4est_connectivity_new_unitsquare();
	p4est = p4est_new_ext(MPI_COMM_WORLD, conn, 0, 0, 0, 0, nullptr, &t);
	p4est_refine(p4est, true, refine_fn, nullptr);
	p4est_partition(p4est, true, nullptr);
	ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
	mesh  = p4est_mesh_new_ext(p4est, ghost, 1, 1, P4EST_CONNECT_FULL);
}
