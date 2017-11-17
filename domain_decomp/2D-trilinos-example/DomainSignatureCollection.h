#ifndef DOMAINSIGNATURECOLLECTION_H
#define DOMAINSIGNATURECOLLECTION_H
#include "Side.h"
#include <array>
#include <bitset>
#include <map>
#include <set>
#include <Teuchos_Comm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <string>
#include <vector>
#include <zoltan_cpp.h>
struct AmgxMap {
    int num_neighbors;
    std::vector<int> neighbors;
    std::vector<int> send_sizes;
    std::vector<const int*> send_maps;
    std::vector<int> recv_sizes;
    std::vector<const int*> recv_maps;
    std::set<Teuchos::ArrayRCP<int>> arrays;
    AmgxMap(){}
    AmgxMap(int num_neighbors){
        this->num_neighbors = num_neighbors;
        neighbors.resize(num_neighbors);
        send_sizes.resize(num_neighbors);
        send_maps.resize(num_neighbors);
        recv_sizes.resize(num_neighbors);
        recv_maps.resize(num_neighbors);
    }
    AmgxMap(const AmgxMap& orig, int n);
};
/**
 * @brief A structure that represents a domain and its relation to other domains.
 */
struct DomainSignature {
	/**
	 * @brief The domain's own id
	 */
	int id = -1;

	int refine_level = 1;

	std::array<int, 8> nbr_id           = {{-1, -1, -1, -1, -1, -1, -1, -1}};
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

	friend bool operator<(const DomainSignature &l, const DomainSignature &r)
	{
		return l.id < r.id;
	}
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
            //do nothing
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
            //do nothing
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
	BlockType      type;
	MatrixBlock(int i, int j, Side main, Side aux, std::bitset<4> neumann, bool zero_patch,
	            BlockType type)
	{
		this->i          = i;
		this->j          = j;
		this->flip_i     = flip_i_table[static_cast<int>(main)][static_cast<int>(aux)];
		this->flip_j     = flip_j_table[static_cast<int>(main)][static_cast<int>(aux)];
		this->right      = main == Side::south || main == Side::west;
		this->neumann    = neumann;
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
using DsMemPtr = int &(DomainSignature::*) (Side);
struct Iface {
	std::set<MatrixBlock> getRowBlocks(int *id, DsMemPtr normal);
	bool            y_axis;
	int             id        = -1;
	int             id_global = -1;
	int             id_local  = -1;
	DomainSignature left, right, extra;
	IfaceType       type;

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
	std::set<MatrixBlock> getGlobalColBlocks();
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
    AmgxMap amgxmap;
	int rank;
    Teuchos::RCP<const Teuchos::Comm<int>> comm;
	int matrix_j_low;
	int matrix_j_high;
	int num_pins;
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
	std::map<int, Iface>           ifaces;
	std::map<int, int>             iface_rev_map;
	std::vector<int> iface_map_vec;
	std::vector<int> iface_off_proc_map_vec;
	std::vector<int> iface_off_proc_vec;
	std::vector<int> iface_off_proc_map_vec_send;
	std::vector<int> iface_off_proc_vec_send;
	std::map<int, int> domain_rev_map;
	std::vector<int> domain_map_vec;

	/**
	 * @brief Default empty constructor.
	 */
	DomainSignatureCollection() = default;

	void enumerateIfaces();
	void determineCoarseness();
	void determineAmrLevel();
	void determineXY();
	void divide();
	/**
	 * @brief Generate a grid of domains.
	 *
	 * @param d_x number of domains in the x direction.
	 * @param d_y number of domains in the y direction.
	 * @param rank the rank of the MPI process.
	 */
	DomainSignatureCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm,int d_x, int d_y, int rank);
	DomainSignatureCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm,int d_x, int d_y, int rank, bool amr);
	DomainSignatureCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm,std::string file_name, int rank);
	/**
	 * @brief Balance the domains over processors using Zoltan
	 */
	void zoltanBalance();
	void zoltanBalanceIfaces();
	void zoltanBalanceDomains();
	/**
	 * @brief Index the interfaces using a Breadth First Search.
	 */
	void indexInterfacesBFS();
	void indexIfacesLocal();
	void indexIfacesGlobal();
	void indexDomainIfacesLocal();
	void setNeumann()
	{
		for (auto &p : domains) {
			p.second.setNeumann();
		}
	}
	void setZeroPatch()
	{
		for (auto &p : domains) {
			p.second.setZeroPatch();
		}
	}
};
struct IfaceZoltanHelper {
	// query functions that respond to requests from Zoltan
	static int get_number_of_objects(void *data, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;

		return dsc->ifaces.size();
	}
	//	static void coord(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR
	// global_id,
	//	                  ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr);
	static void get_object_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID,
	                            ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;

		// In this example, return the IDs of our objects, but no weights.
		// Zoltan will assume equally weighted objects.

		int i = 0;
		for (auto &p : dsc->ifaces) {
			globalID[i]                  = p.first;
			localID[i]                   = p.first;
			std::set<MatrixBlock> blocks = p.second.getGidRowBlocks();
			float                 weight = 1.0;
			for (MatrixBlock b : blocks) {
				if (b.i != b.j) {
					switch (b.type) {
						case BlockType::plain:
							weight += 1.0;
							break;
						case BlockType::fine:
							weight += 1.0;
							break;
						case BlockType::fine_out_left:
							weight += 0.5;
							break;
						case BlockType::fine_out_right:
							weight += 0.5;
							break;
						case BlockType::coarse:
							weight += 1.0;
							break;
						case BlockType::coarse_out_left:
							weight += 1.0;
							break;
						case BlockType::coarse_out_right:
							weight += 1.0;
							break;
					}
				}
			}
			weight      = 1.0;
			obj_wgts[i] = weight;
			i++;
		}
		return;
	}
	static void object_sizes(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *sizes,
	                         int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			sizes[i] = sizeof(dsc->ifaces[global_ids[i]]);
		}
	}
	static void pack_objects(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *dest,
	                         int *sizes, int *idx, char *buf, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			*((Iface *) &buf[idx[i]]) = dsc->ifaces[global_ids[i]];
			dsc->ifaces.erase(global_ids[i]);
		}
	}
	static void unpack_objects(void *data, int num_gid_entries, int num_ids,
	                           ZOLTAN_ID_PTR global_ids, int *sizes, int *idx, char *buf, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			dsc->ifaces[global_ids[i]] = *((Iface *) &buf[idx[i]]);
		}
	}
	static void ZOLTAN_HG_SIZE_CS_FN(void *data, int *num_lists, int *num_pins, int *format,
	                                 int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		*num_lists                     = dsc->ifaces.size();
		*num_pins                      = dsc->num_pins;
		*format                        = ZOLTAN_COMPRESSED_EDGE;
	}
	static void ZOLTAN_HG_CS_FN(void *data, int num_gid_entries, int num_vtx_edge, int num_pins,
	                            int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr,
	                            ZOLTAN_ID_PTR pin_GID, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		int edge_i                     = 0;
		int list_i                     = 0;
		for (auto p : dsc->ifaces) {
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
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;

		return dsc->domains.size();
	}
	//	static void coord(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR
	// global_id,
	//	                  ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr);
	static void get_object_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID,
	                            ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;

		// In this example, return the IDs of our objects, but no weights.
		// Zoltan will assume equally weighted objects.

		int i = 0;
		for (auto &p : dsc->domains) {
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
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			auto ds         = dsc->domains[global_ids[i]];
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
			sizes[i] = sizeof(int) + sizeof(DomainSignature) + num_ifaces * sizeof(Iface);
		}
	}
	static void pack_objects(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *dest,
	                         int *sizes, int *idx, char *buf, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			auto &ds = dsc->domains[global_ids[i]];
			for (int q = 0; q < 8; q++) {
				if (ds.nbr_id[q] != -1) {
					auto &nbr   = dsc->domains[ds.nbr_id[q]];
					nbr.proc[q] = dest[i];
				}
			}
		}
		for (int i = 0; i < num_ids; i++) {
			*((DomainSignature *) &buf[idx[i] + sizeof(int)]) = dsc->domains[global_ids[i]];
			dsc->domains.erase(global_ids[i]);
		}
	}
	static void unpack_objects(void *data, int num_gid_entries, int num_ids,
	                           ZOLTAN_ID_PTR global_ids, int *sizes, int *idx, char *buf, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			dsc->domains[global_ids[i]] = *((DomainSignature *) &buf[idx[i] + sizeof(int)]);
		}
	}
	static int numInterfaces(void *data, int num_gid_entries, int num_lid_entries,
	                         ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
	{
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		auto &ds                       = dsc->domains[*global_id];
		int   num_iface                = 0;
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
		DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
		*ierr                          = ZOLTAN_OK;
		auto &ds                       = dsc->domains[*global_id];
		int   i                        = 0;
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
