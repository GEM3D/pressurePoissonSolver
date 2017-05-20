#ifndef DOMAINSIGNATURECOLLECTION_H
#define DOMAINSIGNATURECOLLECTION_H
#include "Side.h"
#include <array>
#include <bitset>
#include <map>
#include <string>
#include <set>
#include <deque>
/**
 * @brief A structure that represents a domain and its relation to other domains.
 */
struct DomainSignature {
	/**
	 * @brief The domain's own id
	 */
	int id;

	int refine_level = 1;

	std::array<int, 8> nbr_id      = {-1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 8> proc     = {-1, -1, -1, -1, -1, -1, -1, -1};
	std::array<int, 4> global_i = {-1, -1, -1, -1};
	std::array<int, 4> global_i_center = {-1, -1, -1, -1};
	std::array<int, 8> global_i_refined = {-1, -1, -1, -1,-1,-1,-1,-1};
	std::bitset<4> nbr_coarse;
	std::bitset<4> nbr_fine;
	std::bitset<4> left_of_coarse;
	int global_i_west = -1;
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

	friend bool operator<(const DomainSignature &l, const DomainSignature &r)
	{
		return l.id < r.id;
	}
	inline int &index(Side s) { return global_i[static_cast<int>(s)]; }
	inline int &indexCenter(Side s) { return global_i_center[static_cast<int>(s)]; }
	inline int &indexRefinedLeft(Side s) { return global_i_refined[2 * static_cast<int>(s)]; }
	inline int &indexRefinedRight(Side s) { return global_i_refined[1 + 2 * static_cast<int>(s)]; }
	inline int &nbr(Side s) { return nbr_id[2 * static_cast<int>(s)]; }
	inline int &nbrRight(Side s) { return nbr_id[2 * static_cast<int>(s) + 1]; }
	inline bool hasNbr(Side s) const { return nbr_id[static_cast<int>(s) * 2] != -1; }
	inline bool hasFineNbr(Side s) const { return nbr_fine[static_cast<int>(s)]; }
	inline bool hasCoarseNbr(Side s) const { return nbr_coarse[static_cast<int>(s)]; }
	inline bool leftOfCoarse(Side s) const { return left_of_coarse[static_cast<int>(s)]; }
	inline void setHasFineNbr(Side s) { nbr_fine[static_cast<int>(s)] = true; }
	inline void setHasCoarseNbr(Side s) { nbr_coarse[static_cast<int>(s)] = true; }
	inline void setLeftOfCoarse(Side s) { left_of_coarse[static_cast<int>(s)] = true; }
};

/**
 * @brief A collection of Domain Signatures
 *
 * There are two purposes for this class:
 *   -# Partition the domains across processors
 *   -# Determine the number of interfaces, and provide a unique index to each interface.
 */
class DomainSignatureCollection
{
	public:
	int rank;
	/**
	 * @brief Number of total domains.
	 */
	int num_global_domains;
	/**
	 * @brief Number of total interfaces.
	 */
	int num_global_interfaces;
	/**
	 * @brief A map that maps the id of a domain to its domain signature.
	 */
	std::map<int, DomainSignature> domains;

	/**
	 * @brief Default empty constructor.
	 */
	DomainSignatureCollection() = default;

	DomainSignatureCollection(std::string file_name, int rank);
	void determineCoarseness();
	void determineAmrLevel();
	void determineXY();
	/**
	 * @brief Generate a grid of domains.
	 *
	 * @param d_x number of domains in the x direction.
	 * @param d_y number of domains in the y direction.
	 * @param rank the rank of the MPI process.
	 */
	DomainSignatureCollection(int d_x, int d_y, int rank);
	DomainSignatureCollection(int d_x, int d_y, int rank, bool amr);
};
#endif
