#ifndef DOMAINSIGNATURECOLLECTION_H
#define DOMAINSIGNATURECOLLECTION_H
#include "Domain.h"
#include "InterpCase.h"
#include "MyTypeDefs.h"
#include "Side.h"
#include <Teuchos_Comm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <zoltan_cpp.h>
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
	const bool flip_j_table[4][4] = {{false, false, false, false},
	                                 {true, true, true, true},
	                                 {true, true, true, true},
	                                 {false, false, false, false}};
	const bool flip_i_table[4][4] = {{false, false, false, false},
	                                 {true, false, true, false},
	                                 {true, true, true, true},
	                                 {false, true, false, true}};
	int            i, j;
	bool           flip_i, flip_j, right;
	bool           zero_patch;
	std::bitset<4> neumann;
	Side           s;
	InterpCase     type;
	MatrixBlock(int i, int j, Side main, Side aux, std::bitset<4> neumann, bool zero_patch,
	            InterpCase type)
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
	}
	friend bool operator==(const MatrixBlock &l, const MatrixBlock &r)
	{
		return std::tie(l.neumann, l.zero_patch) == std::tie(r.neumann, r.zero_patch);
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
	bool      y_axis;
	int       id        = -1;
	int       id_global = -1;
	int       id_local  = -1;
	Domain    left, right, extra;
	IfaceType type;

	std::vector<int> getEdges()
	{
		std::vector<int> edges;
		for (int i : left.global_i) {
			if (i != -1 && i != id) {
				edges.push_back(i);
			}
		}
		for (int i : right.global_i) {
			if (i != -1 && i != id) {
				edges.push_back(i);
			}
		}
		for (int i : extra.global_i) {
			if (i != -1 && i != id) {
				edges.push_back(i);
			}
		}
		return edges;
	}

	friend bool operator<(const Iface &l, const Iface &r) { return l.id < r.id; }
	std::set<MatrixBlock> getRowBlocks();
	std::set<MatrixBlock> getGlobalRowBlocks();
	std::set<MatrixBlock> getGidRowBlocks();
	std::set<int>         getPins();
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		left.setLocalIndexes(rev_map);
		right.setLocalIndexes(rev_map);
		extra.setLocalIndexes(rev_map);
	}
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		left.setGlobalIndexes(rev_map);
		right.setGlobalIndexes(rev_map);
		extra.setGlobalIndexes(rev_map);
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

	public:
	int                                    rank;
	int                                    n = 4;
	Teuchos::RCP<const Teuchos::Comm<int>> comm;
	int                                    num_pins;
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
	DomainCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm, int d_x, int d_y, int rank);
	DomainCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm, int d_x, int d_y, int rank,
	                 bool amr);
	DomainCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm, std::string file_name, int rank);
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

	void indexDomainIfacesLocal();

	void setNeumann()
	{
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
	Teuchos::RCP<map_type> getSchurRowMap();
	Teuchos::RCP<map_type> getSchurDistMap();
	Teuchos::RCP<map_type> getDomainRowMap();
	int                    getGlobalNumCells() { return num_global_domains * n * n; }
};
struct IfaceZoltanHelper {
	// query functions that respond to requests from Zoltan
	static int get_number_of_objects(void *data, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;

		return dc->ifaces.size();
	}
	//	static void coord(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR
	// global_id,
	//	                  ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr);
	static void get_object_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID,
	                            ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;

		// In this example, return the IDs of our objects, but no weights.
		// Zoltan will assume equally weighted objects.

		int i = 0;
		for (auto &p : dc->ifaces) {
			globalID[i]  = p.first;
			localID[i]   = p.first;
			float weight = 1.0;
			obj_wgts[i]  = weight;
			i++;
		}
		return;
	}
	static void object_sizes(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *sizes,
	                         int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			sizes[i] = sizeof(dc->ifaces[global_ids[i]]);
		}
	}
	static void pack_objects(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *dest,
	                         int *sizes, int *idx, char *buf, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			*((Iface *) &buf[idx[i]]) = dc->ifaces[global_ids[i]];
			dc->ifaces.erase(global_ids[i]);
		}
	}
	static void unpack_objects(void *data, int num_gid_entries, int num_ids,
	                           ZOLTAN_ID_PTR global_ids, int *sizes, int *idx, char *buf, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			dc->ifaces[global_ids[i]] = *((Iface *) &buf[idx[i]]);
		}
	}
	static void ZOLTAN_HG_SIZE_CS_FN(void *data, int *num_lists, int *num_pins, int *format,
	                                 int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		*num_lists           = dc->ifaces.size();
		*num_pins            = dc->num_pins;
		*format              = ZOLTAN_COMPRESSED_EDGE;
	}
	static void ZOLTAN_HG_CS_FN(void *data, int num_gid_entries, int num_vtx_edge, int num_pins,
	                            int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr,
	                            ZOLTAN_ID_PTR pin_GID, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		int edge_i           = 0;
		int list_i           = 0;
		for (auto p : dc->ifaces) {
			std::set<int> list  = p.second.getPins();
			vtxedge_GID[edge_i] = p.first;
			vtxedge_ptr[edge_i] = list_i;
			for (int i : list) {
				pin_GID[list_i] = i;
				list_i++;
			}
			edge_i++;
		}
	}
};

struct DomainZoltanHelper {
	// query functions that respond to requests from Zoltan
	static int get_number_of_objects(void *data, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;

		return dc->domains.size();
	}
	//	static void coord(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR
	// global_id,
	//	                  ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr);
	static void get_object_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID,
	                            ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;

		// In this example, return the IDs of our objects, but no weights.
		// Zoltan will assume equally weighted objects.

		int i = 0;
		for (auto &p : dc->domains) {
			globalID[i] = p.first;
			localID[i]  = p.first;
			i++;
		}
		return;
	}
	static void object_sizes(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *sizes,
	                         int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			auto ds         = dc->domains[global_ids[i]];
			int  num_ifaces = 0;
			if (ds.hasCoarseNbr(Side::north)) {
				num_ifaces += 2;
			} else if (ds.hasFineNbr(Side::north)) {
				num_ifaces += 3;
			} else if (ds.hasNbr(Side::north)) {
				num_ifaces += 1;
			}
			if (ds.hasCoarseNbr(Side::east)) {
				num_ifaces += 2;
			} else if (ds.hasFineNbr(Side::east)) {
				num_ifaces += 3;
			} else if (ds.hasNbr(Side::east)) {
				num_ifaces += 1;
			}
			sizes[i] = sizeof(int) + sizeof(Domain) + num_ifaces * sizeof(Iface);
		}
	}
	static void pack_objects(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *dest,
	                         int *sizes, int *idx, char *buf, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			auto &ds = dc->domains[global_ids[i]];
			for (int q = 0; q < 8; q++) {
				if (ds.nbr_id[q] != -1) {
					auto &nbr   = dc->domains[ds.nbr_id[q]];
					nbr.proc[q] = dest[i];
				}
			}
		}
		for (int i = 0; i < num_ids; i++) {
			*((Domain *) &buf[idx[i] + sizeof(int)]) = dc->domains[global_ids[i]];
			dc->domains.erase(global_ids[i]);
		}
	}
	static void unpack_objects(void *data, int num_gid_entries, int num_ids,
	                           ZOLTAN_ID_PTR global_ids, int *sizes, int *idx, char *buf, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			dc->domains[global_ids[i]] = *((Domain *) &buf[idx[i] + sizeof(int)]);
		}
	}
	static int numInterfaces(void *data, int num_gid_entries, int num_lid_entries,
	                         ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		auto &ds             = dc->domains[*global_id];
		int   num_iface      = 0;
		for (int q = 0; q < 8; q++) {
			if (ds.nbr_id[q] != -1) num_iface++;
		}
		return num_iface;
	}
	static void interfaceList(void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
	                          ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs, int wgt_dim,
	                          float *ewgts, int *ierr)
	{
		DomainCollection *dc = (DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		auto &ds             = dc->domains[*global_id];
		int   i              = 0;
		for (int q = 0; q < 8; q++) {
			if (ds.nbr_id[q] != -1) {
				nbor_global_id[i] = ds.nbr_id[q];
				nbor_procs[i]     = ds.proc[q];
				i++;
			}
		}
	}
	//	static int dimensions(void *data, int *ierr);
};
#endif
