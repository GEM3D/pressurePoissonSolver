#include "DomainSignatureCollection.h"
#include <iostream>
#include <Teuchos_FancyOStream.hpp>
using namespace std;
DomainSignatureCollection::DomainSignatureCollection(int d_x, int d_y, int rank)
{
	num_global_domains    = d_x * d_y;
	num_global_interfaces = 0;
	if (rank == 0) {
		for (int domain_y = 0; domain_y < d_y; domain_y++) {
			for (int domain_x = 0; domain_x < d_x; domain_x++) {
				DomainSignature ds;
				ds.id           = domain_y * d_x + domain_x;
				ds.refine_level = 1;
				if (domain_y != d_y - 1) {
					ds.nbr(Side::north) = (domain_y + 1) * d_x + domain_x;
					ds.proc[0]          = 0;
				}
				if (domain_x != d_x - 1) {
					ds.nbr(Side::east) = domain_y * d_x + domain_x + 1;
					ds.proc[2]         = 0;
				}
				if (domain_y != 0) {
					ds.nbr(Side::south) = (domain_y - 1) * d_x + domain_x;
					ds.proc[4]          = 0;
				}
				if (domain_x != 0) {
					ds.nbr(Side::west) = domain_y * d_x + domain_x - 1;
					ds.proc[6]         = 0;
				}
				ds.x_length    = 1.0 / d_x;
				ds.y_length    = 1.0 / d_y;
				ds.x_start     = 1.0 * domain_x / d_x;
				ds.y_start     = 1.0 * domain_y / d_y;
				domains[ds.id] = ds;
			}
		}
	}
    indexInterfacesBFS();
}
DomainSignatureCollection::DomainSignatureCollection(int d_x, int d_y, int rank,bool amr)
{
	num_global_domains    = d_x * d_y*5;
	num_global_interfaces = 0;
	if (rank == 0) {
		for (int domain_y = 0; domain_y < d_y; domain_y++) {
			for (int domain_x = 0; domain_x < d_x; domain_x++) {
				DomainSignature ds;
				ds.id = domain_y * d_x + domain_x;
                ds.refine_level=1;
				if (domain_y != d_y - 1) {
					ds.nbr(Side::north) = (domain_y + 1) * d_x + domain_x;
					ds.proc[0]          = 0;
				}
				if (domain_x != d_x - 1) {
					ds.nbr(Side::east) = domain_y * d_x + domain_x + 1;
					ds.proc[2]         = 0;
				}
				if (domain_y != 0) {
					ds.nbr(Side::south) = (domain_y - 1) * d_x + domain_x;
					ds.proc[4]          = 0;
				}
				if (domain_x != 0) {
					ds.nbr(Side::west) = domain_y * d_x + domain_x - 1;
					ds.proc[6] = 0;
				}
				ds.x_length    = 1.0 / d_x;
				ds.y_length    = 1.0 / d_y;
				ds.x_start     = 1.0 * domain_x / d_x;
				ds.y_start     = 1.0 * domain_y / d_y;
				domains[ds.id] = ds;
			}
		}
		// create refined grid
		for (int domain_y = 0; domain_y < d_y * 2; domain_y++) {
			for (int domain_x = 0; domain_x < d_x * 2; domain_x++) {
				DomainSignature ds;
				ds.id           = domain_y * d_x * 2 + domain_x + d_x * d_y;
				ds.refine_level = 2;
				if (domain_y != 2 * d_y - 1) {
					ds.nbr(Side::north) = d_x * d_y + 2 * (domain_y + 1) * d_x + domain_x;
					ds.proc[0]          = 0;
				}
				if (domain_x != 2 * d_x - 1) {
					ds.nbr(Side::east) = d_x * d_y + 2 * domain_y * d_x + domain_x + 1;
					ds.proc[2]         = 0;
				}
				if (domain_y != 0) {
					ds.nbr(Side::south) = d_x * d_y + 2 * (domain_y - 1) * d_x + domain_x;
					ds.proc[4]          = 0;
				}
				if (domain_x != 0) {
					ds.nbr(Side::west) = d_x * d_y + 2 * domain_y * d_x + domain_x - 1;
					ds.proc[6]         = 0;
				}
				ds.x_length    = 1.0 / (2 * d_x);
				ds.y_length    = 1.0 / (2 * d_y);
				ds.x_start     = 1.0 + 1.0 * domain_x / (2 * d_x);
				ds.y_start     = 1.0 * domain_y / (2 * d_y);
				domains[ds.id] = ds;
			}
		}
        //stitch together grids
        for(int i=0;i<d_y;i++){
			DomainSignature &left      = domains[i * d_x + d_x - 1];
			DomainSignature &low_nbr   = domains[d_y * d_x + 2 * i * d_x * 2];
			DomainSignature &high_nbr  = domains[d_y * d_x + (2 * i + 1) * d_x * 2];
			left.nbr(Side::east)       = high_nbr.id;
			left.nbrRight(Side::east)  = low_nbr.id;
			left.nbr_fine[1]           = true;
			low_nbr.nbr(Side::west)    = left.id;
			low_nbr.nbr_coarse[3]      = true;
			low_nbr.left_of_coarse[3]  = true;
			high_nbr.nbr(Side::west)   = left.id;
			high_nbr.nbr_coarse[3]     = true;
			high_nbr.left_of_coarse[3] = false;
		}
	}
    indexInterfacesBFS();
}
void DomainSignatureCollection::indexInterfacesBFS()
{
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	int        first = domains.begin()->first;
	queue.push_back(first);
	enqueued.insert(first);
	int curr_i = 0;
	while (!queue.empty()) {
        num_global_interfaces++;
		int              curr = queue.front();
		DomainSignature &d    = domains.at(curr);
		queue.pop_front();
		visited.insert(curr);
        Side s = Side::north;
        do{
			if (d.hasNbr(s) && visited.count(d.nbr(s)) == 0) {
				// a new edge that we have not assigned an index to
				d.index(s) = curr_i;
				curr_i++;

				// fine case
				if (d.hasFineNbr(s)) {
					DomainSignature &nbr_left  = domains.at(d.nbr(s));
					DomainSignature &nbr_right = domains.at(d.nbrRight(s));

					// set center indexes
					nbr_left.indexCenter(!s)  = d.index(s);
					nbr_right.indexCenter(!s) = d.index(s);

					// set left and right indexes index
					nbr_left.index(!s) = curr_i;
					curr_i++;
					nbr_right.index(!s) = curr_i;
					curr_i++;

					// enqueue domains
					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
					if (enqueued.count(d.nbrRight(s)) == 0) {
						queue.push_back(d.nbrRight(s));
						enqueued.insert(d.nbrRight(s));
					}
					// coarse case
				} else if (d.hasCoarseNbr(s)) {
					// TODO
					// normal case
				} else {
					DomainSignature &nbr = domains.at(d.nbr(s));
					nbr.index(!s)        = d.index(s);
					// enqueue domain
					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
				}
			}
			s++;
		} while (s != Side::north);
	}
}
void DomainSignatureCollection::zoltanBalance()
{
	Zoltan *zz = new Zoltan(MPI_COMM_WORLD);

	// parameters
	zz->Set_Param("LB_METHOD", "GRAPH");       /* Zoltan method: "BLOCK" */
	zz->Set_Param("LB_APPROACH", "PARTITION"); /* Zoltan method: "BLOCK" */
	zz->Set_Param("NUM_GID_ENTRIES", "1");     /* global ID is 1 integer */
	zz->Set_Param("NUM_LID_ENTRIES", "1");     /* local ID is 1 integer */
	zz->Set_Param("OBJ_WEIGHT_DIM", "0");      /* we omit object weights */
	zz->Set_Param("AUTO_MIGRATE", "TRUE");     /* we omit object weights */

	// Query functions
	zz->Set_Num_Obj_Fn(DomainSignatureCollection::get_number_of_objects, this);
	zz->Set_Obj_List_Fn(DomainSignatureCollection::get_object_list, this);
	zz->Set_Pack_Obj_Multi_Fn(DomainSignatureCollection::pack_objects, this);
	zz->Set_Unpack_Obj_Multi_Fn(DomainSignatureCollection::unpack_objects, this);
	zz->Set_Obj_Size_Multi_Fn(DomainSignatureCollection::object_sizes, this);
	zz->Set_Num_Edges_Fn(DomainSignatureCollection::numInterfaces, this);
	zz->Set_Edge_List_Fn(DomainSignatureCollection::interfaceList, this);
	// zz->Set_Geom_Fn(DomainSignatureCollection::coord, this);
	// zz->Set_Num_Geom_Fn(DomainSignatureCollection::dimensions, this);

	////////////////////////////////////////////////////////////////
	// Zoltan can now partition the objects in this collection.
	// In this simple example, we assume the number of partitions is
	// equal to the number of processes.  Process rank 0 will own
	// partition 0, process rank 1 will own partition 1, and so on.
	////////////////////////////////////////////////////////////////

	int           changes;
	int           numGidEntries;
	int           numLidEntries;
	int           numImport;
	ZOLTAN_ID_PTR importGlobalIds;
	ZOLTAN_ID_PTR importLocalIds;
	int *         importProcs;
	int *         importToPart;
	int           numExport;
	ZOLTAN_ID_PTR exportGlobalIds;
	ZOLTAN_ID_PTR exportLocalIds;
	int *         exportProcs;
	int *         exportToPart;

	int rc = zz->LB_Partition(changes, numGidEntries, numLidEntries, numImport, importGlobalIds,
	                          importLocalIds, importProcs, importToPart, numExport, exportGlobalIds,
	                          exportLocalIds, exportProcs, exportToPart);

	if (rc != ZOLTAN_OK) {
		cerr << "zoltan error\n";
		delete zz;
		exit(0);
	}
	auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cerr));
	out->setShowProcRank(true);
	*out << "I have " << domains.size() << " domains: ";

    int prev=-100;
    bool range = false;
	for (auto &p : domains)
	{
        int curr = p.second.id;
		if (curr != prev + 1 && !range) {
			*out << curr << "-";
            range = true;
		}else if(curr != prev + 1 && range){
            *out << prev << " " << curr << "-";
        }
		prev = curr;
        
	}

	*out <<prev<< "\n";
}

/*int DomainSignatureCollection::dimensions(void *data, int *ierr)
{
	*ierr = ZOLTAN_OK;
	return 2;
}*/
/*void DomainSignatureCollection::coord(void *data, int num_gid_entries, int num_lid_entries,
                                      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                                      double *geom_vec, int *ierr){
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;

	auto &ds    = dsc->domains[*global_id];
	geom_vec[0] = ds.id % dsc->d_y;
	geom_vec[1] = ds.id % dsc->d_x;
}*/
// query functions that respond to requests from Zoltan
int DomainSignatureCollection::get_number_of_objects(void *data, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;

	return dsc->domains.size();
}

void DomainSignatureCollection::get_object_list(void *data, int sizeGID, int sizeLID,
                                                ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                                int wgt_dim, float *obj_wgts, int *ierr)
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

void DomainSignatureCollection::object_sizes(void *data, int num_gid_entries, int num_lid_entries,
                                             int num_ids, ZOLTAN_ID_PTR global_ids,
                                             ZOLTAN_ID_PTR local_ids, int *sizes, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	for (int i = 0; i < num_ids; i++) {
		sizes[i] = sizeof(dsc->domains[global_ids[i]]);
	}
}
void DomainSignatureCollection::pack_objects(void *data, int num_gid_entries, int num_lid_entries,
                                             int num_ids, ZOLTAN_ID_PTR global_ids,
                                             ZOLTAN_ID_PTR local_ids, int *dest, int *sizes,
                                             int *idx, char *buf, int *ierr)
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
		*((DomainSignature *) &buf[idx[i]]) = dsc->domains[global_ids[i]];
		dsc->domains.erase(global_ids[i]);
	}
}
void DomainSignatureCollection::unpack_objects(void *data, int num_gid_entries, int num_ids,
                                               ZOLTAN_ID_PTR global_ids, int *sizes, int *idx,
                                               char *buf, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	for (int i = 0; i < num_ids; i++) {
		dsc->domains[global_ids[i]] = *((DomainSignature *) &buf[idx[i]]);
	}
}
int DomainSignatureCollection::numInterfaces(void *data, int num_gid_entries, int num_lid_entries,
                                             ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                                             int *ierr)
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
void DomainSignatureCollection::interfaceList(void *data, int num_gid_entries, int num_lid_entries,
                                              ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                                              ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs,
                                              int wgt_dim, float *ewgts, int *ierr)
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
