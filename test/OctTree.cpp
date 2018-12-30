#include "OctTree.h"
#include "catch.hpp"
#include <set>
using namespace std;
TEST_CASE("Tree<3> default constructor works", "[Side]")
{
	Tree<3> tree;
	REQUIRE(tree.nodes.size() == 1);
	REQUIRE(tree.levels.size() == 1);
	REQUIRE(tree.num_levels == 1);
	REQUIRE(tree.levels[1] != nullptr);
	REQUIRE(tree.levels[1]->id == 0);
	REQUIRE(tree.root == 0);
	REQUIRE(tree.max_id == 0);
	// check the root node
	Node<3> n = tree.nodes[0];
	REQUIRE(n.id == 0);
	REQUIRE(n.level == 0);
	REQUIRE(n.parent == -1);
	REQUIRE(n.lengths[0] == 1);
	REQUIRE(n.lengths[1] == 1);
	REQUIRE(n.lengths[2] == 1);
	REQUIRE(n.starts[0] == 0);
	REQUIRE(n.starts[1] == 0);
	REQUIRE(n.starts[2] == 0);
	for (int id : n.nbr_id) {
		REQUIRE(id == -1);
	}
	for (int id : n.child_id) {
		REQUIRE(id == -1);
	}
}
TEST_CASE("Tree<3> refineLeaves() works on single starting node", "[Side]")
{
	Tree<3> tree;
	tree.refineLeaves();
	REQUIRE(tree.nodes.size() == 9);
	REQUIRE(tree.levels.size() == 2);
	REQUIRE(tree.num_levels == 2);
	REQUIRE(tree.levels[1] != nullptr);
	REQUIRE(tree.levels[1]->id >= 0);
	REQUIRE(tree.root == 0);
	REQUIRE(tree.max_id == 8);
	// check the root node
	{
		Node<3> n = tree.nodes[0];
		REQUIRE(n.id == 0);
		REQUIRE(n.level == 0);
		REQUIRE(n.parent == -1);
		REQUIRE(n.lengths[0] == 1);
		REQUIRE(n.lengths[1] == 1);
		REQUIRE(n.lengths[2] == 1);
		REQUIRE(n.starts[0] == 0);
		REQUIRE(n.starts[1] == 0);
		REQUIRE(n.starts[2] == 0);
		set<int> child_ids;
		for (int id : n.nbr_id) {
			REQUIRE(id == -1);
		}
		for (int id : n.child_id) {
			REQUIRE(id != -1);
			child_ids.insert(id);
		}
		REQUIRE(child_ids.size() == 8);
	}
	// check child nodes
	{
		array<int, 8> children = tree.nodes[0].child_id;
		for (Orthant<3> o : Orthant<3>::getValues()) {
			Node<3> child = tree.nodes[children[o.toInt()]];
			REQUIRE(child.id == children[o.toInt()]);
			for (Side<3> s : o.getInteriorSides()) {
				REQUIRE(child.hasNbr(s));
				REQUIRE(child.nbrId(s) == children[o.getInteriorNbrOnSide(s).toInt()]);
			}
			for (Side<3> s : o.getExteriorSides()) {
				REQUIRE(!child.hasNbr(s));
			}
			for (int id : child.child_id) {
				REQUIRE(id == -1);
			}
		}
	}
}
TEST_CASE("Tree<3> refineLeaves() works on single starting node with two calls", "[Side]")
{
	Tree<3> tree;
	tree.refineLeaves();
	tree.refineLeaves();
	REQUIRE(tree.nodes.size() == 73);
	REQUIRE(tree.levels.size() == 3);
	REQUIRE(tree.num_levels == 3);
	REQUIRE(tree.levels.at(1) != nullptr);
	REQUIRE(tree.levels.at(1)->id >= 0);
	REQUIRE(tree.levels.at(2)->id < 9);
	REQUIRE(tree.levels.at(3)->id > 8);
	REQUIRE(tree.root == 0);
	REQUIRE(tree.max_id == 72);
	// check the root node
	{
		Node<3> n = tree.nodes.at(0);
		REQUIRE(n.id == 0);
		REQUIRE(n.level == 0);
		REQUIRE(n.parent == -1);
		REQUIRE(n.lengths[0] == 1);
		REQUIRE(n.lengths[1] == 1);
		REQUIRE(n.lengths[2] == 1);
		REQUIRE(n.starts[0] == 0);
		REQUIRE(n.starts[1] == 0);
		REQUIRE(n.starts[2] == 0);
		set<int> child_ids;
		for (int id : n.nbr_id) {
			REQUIRE(id == -1);
		}
		for (int id : n.child_id) {
			REQUIRE(id != -1);
			child_ids.insert(id);
		}
		REQUIRE(child_ids.size() == 8);
	}
	// check child nodes
	{
		array<int, 8> children = tree.nodes.at(0).child_id;
		for (Orthant<3> o : Orthant<3>::getValues()) {
			Node<3> child = tree.nodes.at(children[o.toInt()]);
			REQUIRE(child.id == children[o.toInt()]);
			for (Side<3> s : o.getInteriorSides()) {
				REQUIRE(child.hasNbr(s));
				REQUIRE(child.nbrId(s) == children[o.getInteriorNbrOnSide(s).toInt()]);
			}
			for (Side<3> s : o.getExteriorSides()) {
				REQUIRE(!child.hasNbr(s));
			}
			// check the child's child nodes
			for (int i = 0; i < 8; i++) {
				Orthant<3> child_o     = i;
				int        id          = child.child_id[i];
				Node<3>    child_child = tree.nodes.at(id);
				REQUIRE(id > 8);
				REQUIRE(id == child_child.id);

				// check interior neighbors
				for (Side<3> s : child_o.getInteriorSides()) {
					REQUIRE(child_child.hasNbr(s));
					REQUIRE(child_child.nbrId(s)
					        == child.child_id[child_o.getInteriorNbrOnSide(s).toInt()]);
				}
				// check exterior neighbors
				{
					auto o_ext       = o.getExteriorSides();
					auto child_o_ext = child_o.getExteriorSides();
					for (int i = 0; i < 3; i++) {
						if (o_ext[i] == child_o_ext[i]) {
							REQUIRE(!child.hasNbr(child_o_ext[i]));
						} else {
							Side<3> s = child_o_ext[i];
							REQUIRE(child_child.hasNbr(s));
							REQUIRE(child_child.nbrId(s)
							        == tree.nodes.at(children[o.getInteriorNbrOnSide(s).toInt()])
							           .child_id[child_o.getExteriorNbrOnSide(s).toInt()]);
						}
					}
				}
				// check children
				for (int id : child_child.child_id) {
					REQUIRE(id == -1);
				}
			}
		}
	}
}
