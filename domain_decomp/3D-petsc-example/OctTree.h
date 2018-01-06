#ifndef OCTREE_H
#define OCTREE_H
#include "OctNode.h"
#include <map>
#include <string>
struct OctTree {
	std::map<int, OctNode> nodes;
	int root;
	OctTree(std::string file_name);
};
#endif
