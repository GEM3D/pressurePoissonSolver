#ifndef OCTNODE_H
#define OCTNODE_H
#include "Side.h"
#include <array>
/**
 * @brief Represents a node in an OctTree
 */
struct OctNode {
	/**
	 * @brief The id of the node.
	 */
	int id = -1;
	/**
	 * @brief The level of the tree that the node lies on.
	 */
	int level = -1;
	/**
	 * @brief The id of the node's parent. -1 if no parent.
	 */
	int parent = -1;
	/**
	 * @brief The x length of this node.
	 */
	double x_length = -1;
	/**
	 * @brief The y length of this node.
	 */
	double y_length = -1;
	/**
	 * @brief The z length of this node.
	 */
	double z_length = -1;
	/**
	 * @brief The x-coordinate for the bottom-south-west of the node.
	 */
	double x_start = -1;
	/**
	 * @brief The y-coordinate for the bottom-south-west of the node.
	 */
	double y_start = -1;
	/**
	 * @brief The z-coordinate for the bottom-south-west of the node.
	 */
	double z_start = -1;
	/**
	 * @brief Array of neighbor ids, -1 if no neighbor
	 */
	std::array<int, 6> nbr_id = {{-1, -1, -1, -1, -1, -1}};
	/**
	 * @brief Array of child ids, -1 if no children
	 */
	std::array<int, 8> child_id = {{-1, -1, -1, -1, -1, -1, -1, -1}};
	/**
	 * @brief Create new node with all values initialized to -1
	 */
	OctNode() = default;
	/**
	 * @brief Create new node that that is the child of a parent and lies on a particular octant of
	 * that parent.
	 *
	 * @param parent The parent of the new node.
	 * @param o the octant of the parent that this node lies on.
	 */
	OctNode(OctNode parent, Octant o)
	{
		this->parent = parent.id;
		level        = parent.level + 1;
		x_length     = parent.x_length / 2;
		y_length     = parent.y_length / 2;
		z_length     = parent.z_length / 2;
		x_start      = o.isOnSide(Side::west) ? parent.x_start : (parent.x_start + this->x_length);
		y_start      = o.isOnSide(Side::south) ? parent.y_start : (parent.y_start + this->y_length);
		z_start = o.isOnSide(Side::bottom) ? parent.z_start : (parent.z_start + this->z_length);
	}
	/**
	 * @return Whether or not this node has any children.
	 */
	bool hasChildren()
	{
		return child_id[0] != -1;
	}
	/**
	 * @brief Return whether or not this node has a neighbor on a particular side.
	 *
	 * @param s the side to check.
	 *
	 * @return The result.
	 */
	bool hasNbr(Side s)
	{
		return nbr_id[s.toInt()] != -1;
	}
	/**
	 * @brief Return reference to the child id for a particular octant.
	 *
	 * @param o the octant
	 *
	 * @return the child id.
	 */
	int &childId(Octant o)
	{
		return child_id[o.toInt()];
	}
	/**
	 * @brief Return reference to neighbor id on a particular side.
	 *
	 * @param s the side
	 *
	 * @return the neighbor id.
	 */
	int &nbrId(Side s)
	{
		return nbr_id[s.toInt()];
	}
};
#endif
