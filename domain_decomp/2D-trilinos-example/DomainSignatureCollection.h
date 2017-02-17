#ifndef DOMAINSIGNATURECOLLECTION_H
#define DOMAINSIGNATURECOLLECTION_H
#include <zoltan_cpp.h>
#include <list>
class DomainSignature
{
    public:
	int id;
	int nbr_north = -1;
	int nbr_south = -1;
	int nbr_east  = -1;
	int nbr_west  = -1;
};

class DomainSignatureCollection
{
	private:
	int  numGlobalObjects=0;
	int  numMyObjects=0;
    std::list<DomainSignature> domains;

	public:
	DomainSignatureCollection() {}
	DomainSignatureCollection(int d_x, int d_y);
	/*	void set_num_global_objects(int n) { numGlobalObjects = n; }
	    int                             get_num_global_objects() { return numGlobalObjects; }
	    void set_num_my_objects(int n) { numMyObjects = n; }
	    int                         get_num_my_objects() { return numMyObjects; }
	    void set_my_global_ids(int *p) { myGlobalIDs = p; }
	    int *                       get_my_global_ids() { return myGlobalIDs; }
	    */
	void zoltanBalance();
	// query functions that respond to requests from Zoltan
	static int get_number_of_objects(void *data, int *ierr);
	static void get_object_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID,
	                            ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr);
	/*static void object_sizes(void *data, int num_gid_entries, int num_lid_entries, int num_ids,
	                         ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *sizes,
	                         int *ierr)
	{
		DommainSigCollection *dsc = (DommainSigCollection *) data;
		*ierr                      = ZOLTAN_OK;
	}*/
};
#endif
