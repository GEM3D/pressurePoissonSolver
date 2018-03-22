#ifndef ZOLTANHELPERS_H
#define ZOLTANHELPERS_H
#include <zoltan.h>
/*
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
*/
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
		*ierr                = ZOLTAN_OK;
		for (int i = 0; i < num_ids; i++) {
			sizes[i] = sizeof(int) + sizeof(Domain);
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
			/*	if (ds.nbr_id[q] != -1) {
					auto &nbr   = dc->domains[ds.nbr_id[q]];
					nbr.proc[q] = dest[i];
				}
                */
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
			//if (ds.nbr_id[q] != -1) num_iface++;
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
		/*	if (ds.nbr_id[q] != -1) {
				nbor_global_id[i] = ds.nbr_id[q];
				nbor_procs[i]     = ds.proc[q];
				i++;
			}*/
		}
	}
	//	static int dimensions(void *data, int *ierr);
};
#endif
