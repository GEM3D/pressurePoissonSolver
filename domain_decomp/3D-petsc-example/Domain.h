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
	int id        = 0;
	int id_local  = 0;
	int id_global = 0;
	int n         = 10;

	int refine_level = 1;

	std::array<int, 6> nbr_id           = {{-1,-1,-1,-1,-1,-1}};
	std::array<int, 6> nbr_id_local     = {{-1,-1,-1,-1,-1,-1}};
	std::array<int, 6> nbr_id_global    = {{-1,-1,-1,-1,-1,-1}};
	std::array<int, 6> proc             = {{-1,-1,-1,-1,-1,-1}};
	std::array<int, 6> g_id             = {{-1,-1,-1,-1,-1,-1}};
	std::array<int, 6> global_i         = {{-1,-1,-1,-1,-1,-1}};
	std::array<int, 6> local_i          = {{-1,-1,-1,-1,-1,-1}};
	std::bitset<6> neumann;
	bool           zero_patch = false;
	/**
	 * @brief The lower left x coordinate of domain
	 */
	double x_start = 0;
	/**
	 * @brief The lower left y coordinate of domain
	 */
	double y_start = 0;
	double z_start = 0;
	/**
	 * @brief length of domain in x direction
	 */
	double x_length = 1;
	/**
	 * @brief length of odmain in y direction
	 */
	double y_length = 1;
	double z_length = 1;

	friend bool operator<(const Domain &l, const Domain &r) { return l.id < r.id; }
	int &gid(Side s) { return g_id[static_cast<int>(s)]; }
	int &globalIndex(Side s) { return global_i[static_cast<int>(s)]; }
	int &index(Side s) { return local_i[static_cast<int>(s)]; }
	inline int &nbr(Side s) { return nbr_id[ static_cast<int>(s)]; }
	inline int &globalNbr(Side s) { return nbr_id_global[2 * static_cast<int>(s)]; }
	inline bool hasNbr(Side s) const { return nbr_id[static_cast<int>(s)] != -1; }
	inline bool isNeumann(Side s) const { return neumann[static_cast<int>(s)]; }
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		for (int i = 0; i < 6; i++) {
			if (g_id[i] != -1) {
				local_i[i] = rev_map.at(g_id[i]);
			}
		}
	}
	void setLocalNeighborIndexes(std::map<int, int> &rev_map)
	{
		id_local = rev_map.at(id);
		for (int i = 0; i < 6; i++) {
			if (nbr_id[i] != -1) {
				nbr_id_local[i] = rev_map.at(nbr_id[i]);
			}
		}
	}

	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		for (int i = 0; i < 6; i++) {
			if (g_id[i] != -1) {
				global_i[i] = rev_map.at(local_i[i]);
			}
		}
	}
	void setGlobalNeighborIndexes(std::map<int, int> &rev_map)
	{
		id_global = rev_map.at(id_local);
		for (int i = 0; i < 6; i++) {
			if (nbr_id_local[i] != -1) {
				nbr_id_global[i] = rev_map.at(nbr_id_local[i]);
			}
		}
	}
	void setNeumann()
	{
		for (int q = 0; q < 6; q++) {
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
