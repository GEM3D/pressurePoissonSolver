#ifndef OCTREE_H
#define OCTREE_H
#include "OctNode.h"
#include <map>
#include <string>
/**
 * @brief Represents an oct-tree.
 */
struct OctTree {
	/**
	 * @brief Map of node id to node objects.
	 */
	std::map<int, OctNode> nodes;
	/**
	 * @brief Map of level to a node object on that level. With 0 begin the root of the tree.
	 */
	std::map<int, OctNode *> levels;
	/**
	 * @brief The id of the root.
	 */
	int root;
	/**
	 * @brief The number of levels in this tree.
	 */
	int num_levels;
	/**
	 * @brief The maximum id value used in this tree.
	 */
	int max_id;
	/**
	 * @brief Create new tree with only a root node.
	 */
	OctTree();
	/**
	 * @brief Read an OctTree from a file.
	 *
	 * @param file_name the file to read from.
	 */
	OctTree(std::string file_name);
	/**
	 * @brief Refine the tree to add one more level.
	 */
	void refineLeaves();
	/**
	 * @brief Refine a particular node.
	 *
	 * @param n the node to refine.
	 */
	void refineNode(OctNode &n);
};
#endif
