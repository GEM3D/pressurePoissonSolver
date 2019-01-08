#ifndef OCTNODE_H
#define OCTNODE_H
#include "Side.h"
#include <array>
/**
 * @brief Represents a node in an OctTree
 */
template <size_t D>
struct Node {
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
	 * @brief The lengths of this node.
	 */
	std::array<double, D> lengths;
	/**
	 * @brief The coordinate for the bottom-south-west of the node.
	 */
	std::array<double, D> starts;
	/**
	 * @brief Array of neighbor ids, -1 if no neighbor
	 */
	std::array<int, Side<D>::num_sides> nbr_id;
	/**
	 * @brief Array of child ids, -1 if no children
	 */
	std::array<int, Orthant<D>::num_orthants> child_id;
	/**
	 * @brief Create new node with all values initialized to -1
	 */
	Node()
	{
		lengths.fill(-1);
		starts.fill(-1);
		nbr_id.fill(-1);
		child_id.fill(-1);
	}
	/**
	 * @brief Create new node that that is the child of a parent and lies on a particular octant of
	 * that parent.
	 *
	 * @param parent The parent of the new node.
	 * @param o the octant of the parent that this node lies on.
	 */
	Node(Node parent, Orthant<D> o)
	{
		this->parent = parent.id;
		level        = parent.level + 1;
		for (size_t i = 0; i < D; i++) {
			lengths[i] = parent.lengths[i] / 2;
			starts[i]
			= o.isOnSide(2 * i) ? parent.starts[i] : (parent.starts[i] + this->lengths[i]);
		}
	}
	/**
	 * @return Whether or not this node has any children.
	 */
	bool hasChildren() const
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
	bool hasNbr(Side<D> s) const
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
	int &childId(Orthant<D> o)
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
	int &nbrId(Side<D> s)
	{
		return nbr_id[s.toInt()];
	}
	int nbrId(Side<D> s) const
	{
		return nbr_id[s.toInt()];
	}
};
#endif
