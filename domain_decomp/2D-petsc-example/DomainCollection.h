#ifndef DOMAINSIGNATURECOLLECTION_H
#define DOMAINSIGNATURECOLLECTION_H
#include "Domain.h"
#include "InterpCase.h"
#include "PW.h"
#include "Side.h"
#include <functional>
#include <map>
#include <memory>
#include <petscvec.h>
#include <set>
#include <string>
#include <vector>
enum class BlockType {
	plain,
	fine,
	fine_out_left,
	fine_out_right,
	coarse,
	coarse_out_left,
	coarse_out_right
};
struct MatrixBlock {
	const bool     flip_j_table[4][4] = {{false, false, false, false},
                                     {true, true, true, true},
                                     {true, true, true, true},
                                     {false, false, false, false}};
	const bool     flip_i_table[4][4] = {{false, false, false, false},
                                     {true, false, true, false},
                                     {true, true, true, true},
                                     {false, true, false, true}};
	int            i, j;
	bool           flip_i, flip_j, right;
	bool           zero_patch;
	std::bitset<4> neumann;
	Side           s;
	InterpCase     type;
	double         length;
	MatrixBlock(int i, int j, Side main, Side aux, std::bitset<4> neumann, bool zero_patch,
	            InterpCase type, double length)
	{
		this->i      = i;
		this->j      = j;
		this->flip_i = flip_i_table[static_cast<int>(main)][static_cast<int>(aux)];
		this->flip_j = flip_j_table[static_cast<int>(main)][static_cast<int>(aux)];
		this->right  = main == Side::south || main == Side::west;
		for (int q = 0; q < 4; q++) {
			this->neumann[q] = neumann[(static_cast<int>(main) + q) % 4];
		}
		this->zero_patch = zero_patch;
		this->type       = type;
		this->s          = aux;
		this->length     = length;
	}
	friend bool operator==(const MatrixBlock &l, const MatrixBlock &r)
	{
		return std::tie(l.neumann, l.zero_patch, l.length)
		       == std::tie(r.neumann, r.zero_patch, r.length);
	}
	friend bool operator<(const MatrixBlock &l, const MatrixBlock &r)
	{
		return std::tie(l.j, l.i, l.right) < std::tie(r.j, r.i, r.right);
	}
};

enum class IfaceType {
	normal,
	coarse_on_left,
	coarse_on_right,
	refined_on_left_left_of_coarse,
	refined_on_left_right_of_coarse,
	refined_on_right_left_of_coarse,
	refined_on_right_right_of_coarse
};
using DsMemPtr = int &(Domain::*) (Side);
struct Iface {
	std::set<MatrixBlock> getRowBlocks(int *id, DsMemPtr normal);
	bool                  y_axis;
	int                   id        = -1;
	int                   id_global = -1;
	int                   id_local  = -1;
	Domain                left, right, extra;
	IfaceType             type;

	std::vector<int> getEdges()
	{
		std::vector<int> edges;
		for (int i : left.global_i) {
			if (i != -1 && i != id) { edges.push_back(i); }
		}
		for (int i : right.global_i) {
			if (i != -1 && i != id) { edges.push_back(i); }
		}
		for (int i : extra.global_i) {
			if (i != -1 && i != id) { edges.push_back(i); }
		}
		return edges;
	}

	friend bool           operator<(const Iface &l, const Iface &r) { return l.id < r.id; }
	std::set<MatrixBlock> getRowBlocks();
	std::set<MatrixBlock> getGlobalRowBlocks();
	std::set<MatrixBlock> getGidRowBlocks();
	std::set<int>         getPins();
	void                  setLocalIndexes(std::map<int, int> &rev_map)
	{
		left.setLocalIndexes(rev_map);
		right.setLocalIndexes(rev_map);
		extra.setLocalIndexes(rev_map);
	}
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		id_global = rev_map.at(id_local);
		left.setGlobalIndexes(rev_map);
		right.setGlobalIndexes(rev_map);
		extra.setGlobalIndexes(rev_map);
	}
	void setLocalNeighborIndexes(std::map<int, int> &rev_map)
	{
		left.setLocalNeighborIndexes(rev_map);
		right.setLocalNeighborIndexes(rev_map);
		if (extra.id != -1) { extra.setLocalNeighborIndexes(rev_map); }
	}
	void setGlobalNeighborIndexes(std::map<int, int> &rev_map)
	{
		left.setGlobalNeighborIndexes(rev_map);
		right.setGlobalNeighborIndexes(rev_map);
		if (extra.id != -1) { extra.setGlobalNeighborIndexes(rev_map); }
	}
	void setZeroPatch()
	{
		left.setZeroPatch();
		right.setZeroPatch();
		extra.setZeroPatch();
	}
	void setNeumann()
	{
		left.setNeumann();
		right.setNeumann();
		extra.setNeumann();
	}
};

/**
 * @brief A collection of Domain Signatures
 *
 * There are two purposes for this class:
 *   -# Partition the domains across processors
 *   -# Determine the number of interfaces, and provide a unique index to each interface.
 */
class DomainCollection
{
	private:
	void enumerateIfaces();
	void determineCoarseness();
	void determineAmrLevel();
	void determineXY();
	void zoltanBalanceIfaces();
	void zoltanBalanceDomains();
	void getRankSize();

	public:
	bool neumann = false;
	int  rank;
	int  size;
	int  n = 4;
	int  num_pins;
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
	std::map<int, Domain> domains;
	std::map<int, Iface>  ifaces;

	std::vector<int> iface_map_vec;
	std::vector<int> iface_off_proc_map_vec;
	std::vector<int> iface_dist_map_vec;
	std::vector<int> domain_map_vec;
	std::vector<int> domain_off_proc_map_vec;

	/**
	 * @brief Default empty constructor.
	 */
	DomainCollection() = default;

	/**
	 * @brief Generate a grid of domains.
	 *
	 * @param d_x number of domains in the x direction.
	 * @param d_y number of domains in the y direction.
	 * @param rank the rank of the MPI process.
	 */
	DomainCollection(int d_x, int d_y);
	DomainCollection(int d_x, int d_y, bool amr);
	DomainCollection(std::string file_name);
	/**
	 * @brief Balance the domains over processors using Zoltan
	 */
	void zoltanBalance();
	void divide();
	/**
	 * @brief Index the interfaces using a Breadth First Search.
	 */
	void indexInterfacesBFS();

	void indexIfacesLocal();
	void indexIfacesGlobal();

	void indexDomainsLocal();
	void indexDomainsGlobal();

	void indexDomainIfacesLocal();

	void setNeumann()
	{
		neumann = true;
		for (auto &p : domains) {
			p.second.setNeumann();
		}
		for (auto &p : ifaces) {
			p.second.setNeumann();
		}
	}
	void setZeroPatch()
	{
		for (auto &p : domains) {
			p.second.setZeroPatch();
		}
		for (auto &p : ifaces) {
			p.second.setZeroPatch();
		}
	}
	PW_explicit<Vec> getNewSchurVec();
	PW_explicit<Vec> getNewSchurDistVec();
	PW_explicit<Vec> getNewDomainVec();

	int    getGlobalNumCells() { return num_global_domains * n * n; }
	double integrate(const Vec u);
	double area();
};
#endif
