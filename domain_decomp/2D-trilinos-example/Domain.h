#ifndef DOMAIN_H
#define DOMAIN_H
#include "Side.h"
#include <array>
#include <bitset>
#include <map>
/**
 * @brief A structure that represents a domain and its relation to other domains.
 */
struct Domain {
	/**
	 * @brief The domain's own id
	 */
	int id        = -1;
	int id_local  = -1;
	int id_global = -1;
	int n         = 10;

	int refine_level = 1;

	std::array<int, 8> nbr_id           = {{-1, -1, -1, -1, -1, -1, -1, -1}};
	std::array<int, 8> nbr_id_local     = {{-1, -1, -1, -1, -1, -1, -1, -1}};
	std::array<int, 8> proc             = {{-1, -1, -1, -1, -1, -1, -1, -1}};
	std::array<int, 4> g_id             = {{-1, -1, -1, -1}};
	std::array<int, 4> g_id_center      = {{-1, -1, -1, -1}};
	std::array<int, 8> g_id_refined     = {{-1, -1, -1, -1, -1, -1, -1, -1}};
	std::array<int, 4> global_i         = {{-1, -1, -1, -1}};
	std::array<int, 4> global_i_center  = {{-1, -1, -1, -1}};
	std::array<int, 8> global_i_refined = {{-1, -1, -1, -1, -1, -1, -1, -1}};
	std::array<int, 4> local_i          = {{-1, -1, -1, -1}};
	std::array<int, 4> local_i_center   = {{-1, -1, -1, -1}};
	std::array<int, 8> local_i_refined  = {{-1, -1, -1, -1, -1, -1, -1, -1}};
	std::bitset<4> nbr_coarse;
	std::bitset<4> nbr_fine;
	std::bitset<4> left_of_coarse;
	std::bitset<4> neumann;
	bool           zero_patch = false;
	/**
	 * @brief The lower left x coordinate of domain
	 */
	double x_start = 0;
	/**
	 * @brief The lower left y coordinate of domain
	 */
	double y_start = 0;
	/**
	 * @brief length of domain in x direction
	 */
	double x_length = 1;
	/**
	 * @brief length of odmain in y direction
	 */
	double y_length = 1;

	friend bool operator<(const Domain &l, const Domain &r) { return l.id < r.id; }
	int &gid(Side s) { return g_id[static_cast<int>(s)]; }
	int &gidCenter(Side s) { return g_id_center[static_cast<int>(s)]; }
	int &gidRefinedLeft(Side s) { return g_id_refined[2 * static_cast<int>(s)]; }
	int &gidRefinedRight(Side s) { return g_id_refined[1 + 2 * static_cast<int>(s)]; }
	int &globalIndex(Side s) { return global_i[static_cast<int>(s)]; }
	int &globalIndexCenter(Side s) { return global_i_center[static_cast<int>(s)]; }
	int &globalIndexRefinedLeft(Side s) { return global_i_refined[2 * static_cast<int>(s)]; }
	int &globalIndexRefinedRight(Side s) { return global_i_refined[1 + 2 * static_cast<int>(s)]; }
	int &index(Side s) { return local_i[static_cast<int>(s)]; }
	int &indexCenter(Side s) { return local_i_center[static_cast<int>(s)]; }
	int &indexRefinedLeft(Side s) { return local_i_refined[2 * static_cast<int>(s)]; }
	int &indexRefinedRight(Side s) { return local_i_refined[1 + 2 * static_cast<int>(s)]; }
	inline int &nbr(Side s) { return nbr_id[2 * static_cast<int>(s)]; }
	inline int &nbrRight(Side s) { return nbr_id[2 * static_cast<int>(s) + 1]; }
	inline bool hasNbr(Side s) const { return nbr_id[static_cast<int>(s) * 2] != -1; }
	inline bool hasFineNbr(Side s) const { return nbr_fine[static_cast<int>(s)]; }
	inline bool hasCoarseNbr(Side s) const { return nbr_coarse[static_cast<int>(s)]; }
	inline bool leftOfCoarse(Side s) const { return left_of_coarse[static_cast<int>(s)]; }
	inline bool isNeumann(Side s) const { return neumann[static_cast<int>(s)]; }
	inline void setHasFineNbr(Side s) { nbr_fine[static_cast<int>(s)] = true; }
	inline void setHasCoarseNbr(Side s) { nbr_coarse[static_cast<int>(s)] = true; }
	inline void setLeftOfCoarse(Side s) { left_of_coarse[static_cast<int>(s)] = true; }
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		for (int i = 0; i < 4; i++) {
			if (g_id[i] != -1) {
				local_i[i] = rev_map.at(g_id[i]);
			}
		}
		try {
			for (int i = 0; i < 4; i++) {
				if (g_id_center[i] != -1) {
					local_i_center[i] = rev_map.at(g_id_center[i]);
				}
			}
			for (int i = 0; i < 8; i++) {
				if (g_id_refined[i] != -1) {
					local_i_refined[i] = rev_map.at(g_id_refined[i]);
				}
			}
		} catch (std::out_of_range oor) {
			// do nothing
		}
	}
	void setLocalNeighborIndexes(std::map<int, int> &rev_map)
	{
		id_local = rev_map.at(id);
		for (int i = 0; i < 8; i++) {
			if (nbr_id[i] != -1) {
				nbr_id_local[i] = rev_map.at(nbr_id[i]);
			}
		}
	}
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		for (int i = 0; i < 4; i++) {
			if (g_id[i] != -1) {
				global_i[i] = rev_map.at(local_i[i]);
			}
		}
		try {
			for (int i = 0; i < 4; i++) {
				if (g_id_center[i] != -1) {
					global_i_center[i] = rev_map.at(local_i_center[i]);
				}
			}
			for (int i = 0; i < 8; i++) {
				if (g_id_refined[i] != -1) {
					global_i_refined[i] = rev_map.at(local_i_refined[i]);
				}
			}
		} catch (std::out_of_range oor) {
			// do nothing
		}
	}
	void setNeumann()
	{
		for (int q = 0; q < 4; q++) {
			neumann[q] = (global_i[q] == -1);
		}
	}
	std::bitset<4> neumannRelative(Side s)
	{
		std::bitset<4> ret;
		for (int q = 0; q < 4; q++) {
			ret[q] = neumann[(static_cast<int>(s) + q) % 4];
		}
		return ret;
	}
	void setZeroPatch()
	{
		if (id == 0) {
			zero_patch = true;
		}
	}
};
#endif
