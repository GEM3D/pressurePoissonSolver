#include "OctTree.h"
#include <fstream>
#include <iostream>
using namespace std;
OctTree::OctTree(string file_name)
{
	ifstream input(file_name, ios_base::binary);
	int      num_nodes;
	int      num_trees;
	input.read((char *) &num_nodes, 4);
	input.read((char *) &num_trees, 4);
	for (int i = 0; i < num_nodes; i++) {
		OctNode n;
		// id level parent
		input.read((char *) &n.id, 3 * 4);
		// length and starts
		input.read((char *) &n.x_length, 6 * 8);

		input.read((char *) &n.nbr_id[0], 6 * 4);
		input.read((char *) &n.child_id[0], 8 * 4);

		if (i == 0) {
			root = n.id;
		}
		nodes[n.id] = n;

		if (n.level > num_levels) {
			num_levels = n.level;
		}
        levels[n.level]=&nodes[n.id];
	}
}
