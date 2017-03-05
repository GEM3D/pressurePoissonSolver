#include "DomainSignatureCollection.h"
#include <iostream>
#include <Teuchos_FancyOStream.hpp>
using namespace std;
DomainSignatureCollection::DomainSignatureCollection(int d_x, int d_y, int rank)
{
	this->d_x = d_x;
	this->d_y = d_y;
	if (rank == 0) {
		for (int domain_y = 0; domain_y < d_y; domain_y++) {
			for (int domain_x = 0; domain_x < d_x; domain_x++) {
				DomainSignature ds;
				ds.id = domain_y * d_x + domain_x;
				if (domain_y != d_y - 1) {
					ds.nbr_north  = (domain_y + 1) * d_x + domain_x;
					ds.proc_north = 0;
				}
				if (domain_x != d_x - 1) {
					ds.nbr_east  = domain_y * d_x + domain_x + 1;
					ds.proc_east = 0;
				}
				if (domain_y != 0) {
					ds.nbr_south  = (domain_y - 1) * d_x + domain_x;
					ds.proc_south = 0;
				}
				if (domain_x != 0) {
					ds.nbr_west  = domain_y * d_x + domain_x - 1;
					ds.proc_west = 0;
				}
				domains[ds.id] = ds;
			}
		}
	}
    indexDomainsBFS();
}
void DomainSignatureCollection::indexDomainsBFS()
{
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	int        first = domains.begin()->first;
	queue.push_back(first);
	enqueued.insert(first);
	int curr_i = 0;
	while (!queue.empty()) {
		int              curr = queue.front();
		DomainSignature &d    = domains.at(curr);
		queue.pop_front();
		visited.insert(curr);
		if (d.nbr_north != -1 && visited.count(d.nbr_north) == 0) {
			// a new edge that we have not assigned an index to
			DomainSignature &nbr = domains.at(d.nbr_north);
			d.global_i_north     = curr_i;
			nbr.global_i_south   = curr_i;
			if (enqueued.count(d.nbr_north) == 0) {
				queue.push_back(d.nbr_north);
				enqueued.insert(d.nbr_north);
			}
			curr_i++;
		}
		if (d.nbr_east != -1 && visited.count(d.nbr_east) == 0) {
			// a new edge that we have not assigned an index to
			DomainSignature &nbr = domains.at(d.nbr_east);
			d.global_i_east      = curr_i;
			nbr.global_i_west    = curr_i;
			if (enqueued.count(d.nbr_east) == 0) {
				queue.push_back(d.nbr_east);
				enqueued.insert(d.nbr_east);
			}
			curr_i++;
		}
		if (d.nbr_south != -1 && visited.count(d.nbr_south) == 0) {
			// a new edge that we have not assigned an index to
			DomainSignature &nbr = domains.at(d.nbr_south);
			d.global_i_south     = curr_i;
			nbr.global_i_north   = curr_i;
			if (enqueued.count(d.nbr_south) == 0) {
				queue.push_back(d.nbr_south);
				enqueued.insert(d.nbr_south);
			}
			curr_i++;
		}
		if (d.nbr_west != -1 && visited.count(d.nbr_west) == 0) {
			// a new edge that we have not assigned an index to
			DomainSignature &nbr = domains.at(d.nbr_west);
			d.global_i_west      = curr_i;
			nbr.global_i_east    = curr_i;
			if (enqueued.count(d.nbr_west) == 0) {
				queue.push_back(d.nbr_west);
				enqueued.insert(d.nbr_west);
			}
			curr_i++;
		}
	}
}
void DomainSignatureCollection::zoltanBalance()
{
	Zoltan *zz = new Zoltan(MPI_COMM_WORLD);

	// parameters
	zz->Set_Param("LB_METHOD", "GRAPH");   /* Zoltan method: "BLOCK" */
	zz->Set_Param("LB_APPROACH", "PARTITION");   /* Zoltan method: "BLOCK" */
	zz->Set_Param("NUM_GID_ENTRIES", "1"); /* global ID is 1 integer */
	zz->Set_Param("NUM_LID_ENTRIES", "1"); /* local ID is 1 integer */
	zz->Set_Param("OBJ_WEIGHT_DIM", "0");  /* we omit object weights */
	zz->Set_Param("AUTO_MIGRATE", "TRUE");  /* we omit object weights */

	// Query functions
	zz->Set_Num_Obj_Fn(DomainSignatureCollection::get_number_of_objects, this);
	zz->Set_Obj_List_Fn(DomainSignatureCollection::get_object_list, this);
	zz->Set_Pack_Obj_Multi_Fn(DomainSignatureCollection::pack_objects, this);
	zz->Set_Unpack_Obj_Multi_Fn(DomainSignatureCollection::unpack_objects, this);
	zz->Set_Obj_Size_Multi_Fn(DomainSignatureCollection::object_sizes, this);
	zz->Set_Num_Edges_Fn(DomainSignatureCollection::numInterfaces, this);
	zz->Set_Edge_List_Fn(DomainSignatureCollection::interfaceList, this);
	zz->Set_Geom_Fn(DomainSignatureCollection::coord, this);
	zz->Set_Num_Geom_Fn(DomainSignatureCollection::dimensions, this);

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

int DomainSignatureCollection::dimensions(void *data, int *ierr)
{
	*ierr = ZOLTAN_OK;
	return 2;
}
void DomainSignatureCollection::coord(void *data, int num_gid_entries, int num_lid_entries,
                                      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                                      double *geom_vec, int *ierr){
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;

	auto &ds    = dsc->domains[*global_id];
	geom_vec[0] = ds.id % dsc->d_y;
	geom_vec[1] = ds.id % dsc->d_x;
}
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
		if (wgt_dim == 1) {
			int weight = 0;
			if (p.second.nbr_north != -1) weight += 1;
			if (p.second.nbr_east != -1) weight += 1;
			obj_wgts[i] = weight;
		}
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
		if (ds.nbr_north != -1) {
			auto &nbr          = dsc->domains[ds.nbr_north];
			nbr.proc_south= dest[i];
		}
		if (ds.nbr_east != -1) {
			auto &nbr         = dsc->domains[ds.nbr_east];
			nbr.proc_west= dest[i];
		}
		if (ds.nbr_south != -1) {
			auto &nbr          = dsc->domains[ds.nbr_south];
			nbr.proc_north= dest[i];
		}
		if (ds.nbr_west != -1) {
			auto &nbr         = dsc->domains[ds.nbr_west];
			nbr.proc_east= dest[i];
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
                                              int *ierr){
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	auto &ds                       = dsc->domains[*global_id];
    int num_iface=0;
	if (ds.nbr_north != -1) num_iface++;
	if (ds.nbr_east != -1) num_iface++;
	if (ds.nbr_south != -1) num_iface++;
	if (ds.nbr_west != -1) num_iface++;
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
    int i=0;
	if (ds.nbr_north != -1){
        nbor_global_id[i] = ds.nbr_north;
        nbor_procs[i] = ds.proc_north;
        i++;
    }
	if (ds.nbr_east != -1){
        nbor_global_id[i] = ds.nbr_east;
        nbor_procs[i] = ds.proc_east;
        i++;
    }
	if (ds.nbr_south != -1){
        nbor_global_id[i] = ds.nbr_south;
        nbor_procs[i] = ds.proc_south;
        i++;
    }
	if (ds.nbr_west != -1){
        nbor_global_id[i] = ds.nbr_west;
        nbor_procs[i] = ds.proc_west;
        i++;
    }
}
