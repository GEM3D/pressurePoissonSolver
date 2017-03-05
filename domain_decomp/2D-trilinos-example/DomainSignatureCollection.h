#ifndef DOMAINSIGNATURECOLLECTION_H
#define DOMAINSIGNATURECOLLECTION_H
#include <map>
#include <zoltan_cpp.h>
/**
 * @brief A structure that represents a domain and its relation to other domains.
 */
struct DomainSignature {
	/**
	 * @brief The domain's own id
	 */
	int id;
	/**
	 * @brief The id of the neighbor to it's north.
	 * It should be set to -1 if there is no neighbor.
	 */
	int nbr_north = -1;
	/**
	 * @brief The id of the neighbor to it's south.
	 * It should be set to -1 if there is no neighbor.
	 */
	int nbr_south = -1;
	/**
	 * @brief The id of the neighbor to it's east.
	 * It should be set to -1 if there is no neighbor.
	 */
	int nbr_east = -1;
	/**
	 * @brief The id of the neighbor to it's west.
	 * It should be set to -1 if there is no neighbor.
	 */
	int nbr_west = -1;
	/**
	 * @brief The processor id that it's neighbor to the north resides on.
	 * It should be set to -1 if there is no neighbor.
	 */
	int proc_north = -1;
	/**
	 * @brief The processor id that it's neighbor to the south resides on.
	 * It should be set to -1 if there is no neighbor.
	 */
	int proc_south = -1;
	/**
	 * @brief The processor id that it's neighbor to the east resides on.
	 * It should be set to -1 if there is no neighbor.
	 */
	int proc_east = -1;
	/**
	 * @brief The processor id that it's neighbor to the west resides on.
	 * It should be set to -1 if there is no neighbor.
	 */
	int proc_west = -1;
	/**
	 * @brief The northern interface's index in the interface array.
	 * It should be set to -1 if there is no neighbor.
	 */
	int global_i_north = -1;
	/**
	 * @brief The southern interface's index in the interface array.
	 * It should be set to -1 if there is no neighbor.
	 */
	int global_i_south = -1;
	/**
	 * @brief The eastern interface's index in the interface array.
	 * It should be set to -1 if there is no neighbor.
	 */
	int global_i_east = -1;
	/**
	 * @brief The western interface's index in the interface array.
	 * It should be set to -1 if there is no neighbor.
	 */
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
	double x_length = 0;
	/**
	 * @brief length of odmain in y direction
	 */
	double y_length = 0;

	friend bool operator<(const DomainSignature &l, const DomainSignature &r)
	{
		return l.id < r.id;
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
	/**
	 * @brief Generate a grid of domains.
	 *
	 * @param d_x number of domains in the x direction.
	 * @param d_y number of domains in the y direction.
	 * @param rank the rank of the MPI process.
	 */
	DomainSignatureCollection(int d_x, int d_y, int rank);
    /**
     * @brief Balance the domains over processors using Zoltan
     */
	void zoltanBalance();
    /**
     * @brief Index the interfaces using a Breadth First Search.
     */
	void indexInterfacesBFS();

	// query functions that respond to requests from Zoltan
	static int get_number_of_objects(void *data, int *ierr);
	//	static void coord(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR
	//global_id,
	//	                  ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr);
	static void get_object_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID,
	                            ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr);
	static void object_sizes(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *sizes,
	                         int *ierr);
	static void pack_objects(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *dest,
	                         int *sizes, int *idx, char *buf, int *ierr);
	static void unpack_objects(void *data, int num_gid_entries, int num_ids,
	                           ZOLTAN_ID_PTR global_ids, int *sizes, int *idx, char *buf,
	                           int *ierr);
	static int numInterfaces(void *data, int num_gid_entries, int num_lid_entries,
	                         ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr);
	static void interfaceList(void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
	                          ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs, int wgt_dim,
	                          float *ewgts, int *ierr);
	//	static int dimensions(void *data, int *ierr);
};
#endif
