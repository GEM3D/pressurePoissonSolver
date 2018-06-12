#ifndef OCTNODE_H
#define OCTNODE_H
#include "Side.h"
#include <array>
struct OctNode {
	int                id       = -1;
	int                level    = -1;
	int                parent   = -1;
	double             x_length = -1;
	double             y_length = -1;
	double             z_length = -1;
	double             x_start  = -1;
	double             y_start  = -1;
	double             z_start  = -1;
	std::array<int, 6> nbr_id   = {{-1, -1, -1, -1, -1, -1}};
	std::array<int, 8> child_id = {{-1, -1, -1, -1, -1, -1, -1, -1}};
	bool               hasChildren() { return child_id[0] != -1; }
	bool               hasNbr(Side s) { return nbr_id[static_cast<int>(s)] != -1; }
	int &              nbr(Side s) { return nbr_id[static_cast<int>(s)]; }
	int                child(Oct o) { return child_id[static_cast<int>(o)]; }
};
#endif
