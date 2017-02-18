#include "DomainSignatureCollection.h"
#include <iostream>
#include <Teuchos_FancyOStream.hpp>
using namespace std;
DomainSignatureCollection::DomainSignatureCollection(int d_x, int d_y){
	for (int domain_y = 0; domain_y < d_y; domain_y++) {
		for (int domain_x = 0; domain_x < d_x; domain_x++) {
			DomainSignature ds;
			ds.id = domain_y * d_x + domain_x;
			if (domain_y != d_y - 1) {
				ds.nbr_north = (domain_y + 1) * d_x + domain_x;
			}
			if (domain_x != d_x - 1) {
				ds.nbr_east = domain_y * d_x + domain_x + 1;
			}
			if (domain_y != 0) {
				ds.nbr_north = (domain_y - 1) * d_x + domain_x;
			}
			if (domain_x != d_x - 1) {
				ds.nbr_east = domain_y * d_x + domain_x - 1;
			}
			domains[ds.id] = ds;
		}
	}
}
void DomainSignatureCollection::zoltanBalance()
{
	Zoltan *zz = new Zoltan(MPI_COMM_WORLD);

	// parameters
	zz->Set_Param("LB_METHOD", "BLOCK");   /* Zoltan method: "BLOCK" */
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
	auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
	out->setShowProcRank(true);
	*out << "I have " << domains.size() << " domains: ";

	for (auto &p : domains)
	{
		*out << p.first << " ";
	}

	*out << "\n";
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
                                             int *idx, char *buf, int *ierr){
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	for (int i = 0; i < num_ids; i++) {
		*((DomainSignature*) &buf[idx[i]]) = dsc->domains[global_ids[i]];
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
