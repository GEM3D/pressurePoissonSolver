#ifndef OCTREE_H
#define OCTREE_H
#include "OctNode.h"
#include <map>
#include <string>
struct OctTree {
	std::map<int, OctNode> nodes;
	std::map<int, OctNode*> levels;
	int root;
    int num_levels=1;
    OctTree()=default;
	OctTree(std::string file_name);
};
#endif
