#include "OctTree.h"
#include <deque>
#include <fstream>
#include <iostream>
#include <set>
using namespace std;
OctTree::OctTree()
{
	OctNode new_root;
	new_root.id       = 0;
	new_root.x_length = 1;
	new_root.y_length = 1;
	new_root.z_length = 1;
	new_root.x_start  = 0;
	new_root.y_start  = 0;
	new_root.z_start  = 0;
	new_root.level    = 0;
	nodes[0]          = new_root;
	levels[0]         = &nodes[0];
	root              = 0;
	max_id            = 0;
	num_levels        = 1;
}
OctTree::OctTree(string file_name)
{
	ifstream input(file_name, ios_base::binary);
	int      num_nodes;
	int      num_trees;
	input.read((char *) &num_nodes, 4);
	input.read((char *) &num_trees, 4);
	num_levels = 0;
	max_id = 0;
	for (int i = 0; i < num_nodes; i++) {
		OctNode n;
		// id level parent
		input.read((char *) &n.id, 3 * 4);
		// length and starts
		input.read((char *) &n.x_length, 6 * 8);

		input.read((char *) &n.nbr_id[0], 6 * 4);
		input.read((char *) &n.child_id[0], 8 * 4);

		if (i == 0) { root = n.id; }
		max_id      = max(max_id, n.id);
		nodes[n.id] = n;

		if (n.level > num_levels) { num_levels = n.level; }
		levels[n.level] = &nodes[n.id];
	}
}
void OctTree::refineLeaves()
{
	OctNode root_node = nodes[root];
	OctNode child     = root_node;
	int     level     = 0;
	while (child.hasChildren()) {
		child = nodes[child.child_id[0]];
		level++;
	}
	deque<pair<int, int>> q;
	set<pair<int, int>>   qed;
	q.push_back(make_pair(level, child.id));
	qed.insert(make_pair(level, child.id));

	while (!q.empty()) {
		pair<int, int> p     = q.front();
		OctNode        n     = nodes[p.second];
		int            level = p.first;
		q.pop_front();

		// set and enqueue nbrs
		for (Side s : Side::getValues()) {
			// coarser
			if (n.nbrId(s) == -1 && n.parent != -1 && nodes[n.parent].nbrId(s) != -1) {
				OctNode        parent = nodes[n.parent];
				OctNode        nbr    = nodes[parent.nbrId(s)];
				pair<int, int> p(level-1, nbr.id);
                level--;
				if (!qed.count(p)) {
					q.push_back(p);
					qed.insert(p);
				}
				// finer
			} else if (n.nbrId(s) != -1 && nodes[n.nbrId(s)].hasChildren()) {
				OctNode nbr  = nodes[n.nbrId(s)];
				auto    octs = Octant::getValuesOnSide(s.opposite());
                level++;
				for (int i = 0; i < 4; i++) {
					int            id = nbr.childId(octs[i]);
					pair<int, int> p(level, id);
					if (!qed.count(p)) {
						q.push_back(p);
						qed.insert(p);
					}
				}
				// normal
			} else if (n.nbrId(s) != -1) {
				int            id = n.nbrId(s);
				pair<int, int> p(level, id);
				if (!qed.count(p)) {
					q.push_back(p);
					qed.insert(p);
				}
			}
		}
	}
	for (auto p : qed) {
		refineNode(nodes[p.second]);
	}
	levels[num_levels+1] = &nodes[levels[num_levels]->child_id[0]];
	this->num_levels++;
}
void OctTree::refineNode(OctNode &n)
{
	array<OctNode, 8> new_children;
	for (int i = 0; i < 8; i++) {
		new_children[i] = OctNode(n, i);
		OctNode &child  = new_children[i];
		max_id++;
		child.id      = max_id;
		n.child_id[i] = max_id;
	}
	// set new neighbors
	for (Octant o : Octant::getValues()) {
		for (Side s : o.getInteriorSides()) {
			new_children[o.toInt()].nbrId(s) = new_children[o.getInteriorNbrOnSide(s).toInt()].id;
		}
	}

	// set outer neighbors
	for (Side s : Side::getValues()) {
		if (n.hasNbr(s) && nodes[n.nbrId(s)].hasChildren()) {
			OctNode &nbr = nodes[n.nbrId(s)];
			for (Octant o : Octant::getValuesOnSide(s)) {
				OctNode &child                = new_children[o.toInt()];
				OctNode &nbr_child            = nodes[nbr.childId(o.getExteriorNbrOnSide(s))];
				child.nbrId(s)                = nbr_child.id;
				nbr_child.nbrId(s.opposite()) = child.id;
			}
		}
	}
	// add nodes
	for (OctNode child : new_children) {
		nodes[child.id] = child;
	}
}
