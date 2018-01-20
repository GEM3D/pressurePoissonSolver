#include "DomainCollection.h"
#include "ZoltanHelpers.h"
#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <numeric>
#include <petscao.h>
#include <set>
#include <sstream>
#include <utility>
#include <zoltan.h>
using namespace std;
DomainCollection::DomainCollection(OctTree t, int level)
{
	OctNode root  = t.nodes[t.root];
	OctNode child = root;
	for (int i = 1; i < level; i++) {
		child = t.nodes[child.child_id[0]];
	}
	deque<int> q;
	set<int>   qed;
	q.push_back(child.id);
	qed.insert(child.id);

	while (!q.empty()) {
		Domain  d;
		OctNode n = t.nodes[q.front()];
		q.pop_front();

		d.id        = n.id;
		d.x_length  = n.x_length;
		d.y_length  = n.y_length;
		d.z_length  = n.z_length;
		d.x_start   = n.x_start;
		d.y_start   = n.y_start;
		d.z_start   = n.z_start;
		d.child_id  = n.child_id;
		d.parent_id = n.parent;

		Side s = Side::west;
		// set and enqueue nbrs
		do {
			if (n.nbr(s) != -1) {
				int id = n.nbr(s);
				if (!qed.count(id)) {
					q.push_back(id);
					qed.insert(id);
				}
				d.getNbrInfoPtr(s) = new NormalNbrInfo(id);
			}
			s++;
		} while (s != Side::west);
		domains[d.id] = d;
	}
	num_global_domains = domains.size();
	reIndex();
}
DomainCollection::DomainCollection(OctTree t)
{
	OctNode root  = t.nodes[t.root];
	OctNode child = root;
	while (child.hasChildren()) {
		child = t.nodes[child.child_id[0]];
	}
	deque<int> q;
	set<int>   qed;
	q.push_back(child.id);
	qed.insert(child.id);

	while (!q.empty()) {
		Domain  d;
		OctNode n = t.nodes[q.front()];
		q.pop_front();

		d.id       = n.id;
		d.x_length = n.x_length;
		d.y_length = n.y_length;
		d.z_length = n.z_length;
		d.x_start  = n.x_start;
		d.y_start  = n.y_start;
		d.z_start  = n.z_start;

		Side s = Side::west;
		// set and enqueue nbrs
		do {
			if (n.nbr(s) == -1 && n.parent != -1 && t.nodes[n.parent].nbr(s) != -1) {
				OctNode parent = t.nodes[n.parent];
				OctNode nbr    = t.nodes[parent.nbr(s)];
				auto    octs   = getOctsOnSide(s);
				int     quad   = 0;
				while (parent.child(octs[quad]) != n.id) {
					quad++;
				}
				d.getNbrInfoPtr(s) = new CoarseNbrInfo(nbr.id, quad);
				if (!qed.count(nbr.id)) {
					q.push_back(nbr.id);
					qed.insert(nbr.id);
				}
			} else if (n.nbr(s) != -1 && t.nodes[n.nbr(s)].hasChildren()) {
				OctNode nbr  = t.nodes[n.nbr(s)];
				auto    octs = getOctsOnSide(~s);
				array<int, 4> nbr_ids;
				for (int i = 0; i < 4; i++) {
					int id     = nbr.child(octs[i]);
					nbr_ids[i] = id;
					if (!qed.count(id)) {
						q.push_back(id);
						qed.insert(id);
					}
				}
				d.getNbrInfoPtr(s) = new FineNbrInfo(nbr_ids);
			} else if (n.nbr(s) != -1) {
				int id = n.nbr(s);
				if (!qed.count(id)) {
					q.push_back(id);
					qed.insert(id);
				}
				d.getNbrInfoPtr(s) = new NormalNbrInfo(id);
			}
			s++;
		} while (s != Side::west);
		domains[d.id] = d;
	}
	num_global_domains = domains.size();
	reIndex();
}
DomainCollection::DomainCollection(int d_x, int d_y, int d_z)
{
	auto getID = [&](int x, int y, int z) { return x + y * d_y + z * d_z * d_z; };
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	num_global_domains = d_x * d_y * d_z;
	if (rank == 0) {
		for (int domain_z = 0; domain_z < d_z; domain_z++) {
			for (int domain_y = 0; domain_y < d_y; domain_y++) {
				for (int domain_x = 0; domain_x < d_x; domain_x++) {
					Domain ds;
					ds.id           = getID(domain_x, domain_y, domain_z);
					ds.refine_level = 1;
					if (domain_x != 0) {
						ds.getNbrInfoPtr(Side::west)
						= new NormalNbrInfo(getID(domain_x - 1, domain_y, domain_z));
					}
					if (domain_x != d_x - 1) {
						ds.getNbrInfoPtr(Side::east)
						= new NormalNbrInfo(getID(domain_x + 1, domain_y, domain_z));
					}
					if (domain_y != 0) {
						ds.getNbrInfoPtr(Side::south)
						= new NormalNbrInfo(getID(domain_x, domain_y - 1, domain_z));
					}
					if (domain_y != d_y - 1) {
						ds.getNbrInfoPtr(Side::north)
						= new NormalNbrInfo(getID(domain_x, domain_y + 1, domain_z));
					}
					if (domain_z != 0) {
						ds.getNbrInfoPtr(Side::bottom)
						= new NormalNbrInfo(getID(domain_x, domain_y, domain_z - 1));
					}
					if (domain_z != d_z - 1) {
						ds.getNbrInfoPtr(Side::top)
						= new NormalNbrInfo(getID(domain_x, domain_y, domain_z + 1));
					}
					ds.x_length    = 1.0 / d_x;
					ds.y_length    = 1.0 / d_y;
					ds.z_length    = 1.0 / d_z;
					ds.x_start     = 1.0 * domain_x / d_x;
					ds.y_start     = 1.0 * domain_y / d_y;
					ds.z_start     = 1.0 * domain_z / d_y;
					domains[ds.id] = ds;
				}
			}
		}
	}
	reIndex();
}
void DomainCollection::reIndex() { indexDomainsLocal(); }
void DomainCollection::divide() {}
void DomainCollection::zoltanBalance() { zoltanBalanceDomains(); }
void DomainCollection::zoltanBalanceDomains()
{
	struct Zoltan_Struct *zz = Zoltan_Create(MPI_COMM_WORLD);

	// parameters
	Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");       /* Zoltan method: "BLOCK" */
	Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); /* Zoltan method: "BLOCK" */
	Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");     /* global ID is 1 integer */
	Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");     /* local ID is 1 integer */
	Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");      /* we omit object weights */
	Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");     /* we omit object weights */

	// Query functions
	Zoltan_Set_Num_Obj_Fn(zz, DomainZoltanHelper::get_number_of_objects, this);
	Zoltan_Set_Obj_List_Fn(zz, DomainZoltanHelper::get_object_list, this);
	Zoltan_Set_Pack_Obj_Multi_Fn(zz, DomainZoltanHelper::pack_objects, this);
	Zoltan_Set_Unpack_Obj_Multi_Fn(zz, DomainZoltanHelper::unpack_objects, this);
	Zoltan_Set_Obj_Size_Multi_Fn(zz, DomainZoltanHelper::object_sizes, this);
	Zoltan_Set_Num_Edges_Fn(zz, DomainZoltanHelper::numInterfaces, this);
	Zoltan_Set_Edge_List_Fn(zz, DomainZoltanHelper::interfaceList, this);
	// zz->Set_Geom_Fn(DomainCollection::coord, this);
	// zz->Set_Num_Geom_Fn(DomainCollection::dimensions, this);

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

	int rc = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries, &numImport,
	                             &importGlobalIds, &importLocalIds, &importProcs, &importToPart,
	                             &numExport, &exportGlobalIds, &exportLocalIds, &exportProcs,
	                             &exportToPart);

	if (rc != ZOLTAN_OK) {
		cerr << "zoltan error\n";
		Zoltan_Destroy(&zz);
		exit(0);
	}
	Zoltan_Destroy(&zz);
	cout << "I have " << domains.size() << " domains: ";

#if DD_DEBUG
	int  prev  = -100;
	bool range = false;
	for (auto &p : domains) {
		int curr = p.second.id;
		if (curr != prev + 1 && !range) {
			cout << curr << "-";
			range = true;
		} else if (curr != prev + 1 && range) {
			cout << prev << " " << curr << "-";
		}
		prev = curr;
	}

	cout << prev << "\n";
#endif
	cout << endl;
}
void DomainCollection::indexDomainsGlobal()
{
	// global indices are going to be sequentially increasing with rank
	int local_size = domains.size();
	int start_i;
	MPI_Scan(&local_size, &start_i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	start_i -= local_size;
	vector<int> new_global(local_size);
	iota(new_global.begin(), new_global.end(), start_i);

	// create map for gids
	PW<AO> ao;
	AOCreateMapping(MPI_COMM_WORLD, local_size, &domain_map_vec[0], &new_global[0], &ao);

	// get global indices that we want to recieve for dest vector
	vector<int> inds = domain_map_vec;
	for (int i : domain_off_proc_map_vec) {
		inds.push_back(i);
	}

	// get new global indices
	AOApplicationToPetsc(ao, inds.size(), &inds[0]);
	map<int, int> rev_map;
	for (size_t i = 0; i < inds.size(); i++) {
		rev_map[i] = inds[i];
	}

	for (auto &p : domains) {
		p.second.setGlobalNeighborIndexes(rev_map);
	}
	for (size_t i = 0; i < domain_map_vec.size(); i++) {
		domain_map_vec[i] = inds[i];
	}
	for (size_t i = 0; i < domain_off_proc_map_vec.size(); i++) {
		domain_off_proc_map_vec[i] = inds[domain_map_vec.size() + i];
	}
}
void DomainCollection::indexDomainsLocal()
{
	int         curr_i = 0;
	vector<int> map_vec;
	vector<int> off_proc_map_vec;
	map<int, int> rev_map;
	set<int> offs;
	if (!domains.empty()) {
		set<int> todo;
		for (auto &p : domains) {
			todo.insert(p.first);
		}
		set<int> enqueued;
		while (!todo.empty()) {
			deque<int> queue;
			queue.push_back(*todo.begin());
			enqueued.insert(*todo.begin());
			while (!queue.empty()) {
				int i = queue.front();
				todo.erase(i);
				queue.pop_front();
				map_vec.push_back(i);
				Domain &d  = domains[i];
				rev_map[i] = curr_i;
				d.id_local = curr_i;
				curr_i++;
				for (int i : d.getNbrIds()) {
					if (!enqueued.count(i)) {
						enqueued.insert(i);
						if (domains.count(i)) {
							queue.push_back(i);
						} else {
							if (!offs.count(i)) {
								offs.insert(i);
								off_proc_map_vec.push_back(i);
							}
						}
					}
				}
			}
		}
	}
	// map off proc
	for (int i : off_proc_map_vec) {
		rev_map[i] = curr_i;
		curr_i++;
	}
	for (auto &p : domains) {
		p.second.setLocalNeighborIndexes(rev_map);
	}
	// domain_rev_map          = rev_map;
	domain_map_vec          = map_vec;
	domain_off_proc_map_vec = off_proc_map_vec;
	indexDomainsGlobal();
}
double DomainCollection::integrate(const Vec u)
{
	double  sum = 0;
	double *u_view;
	VecGetArray(u, &u_view);

	for (auto &p : domains) {
		Domain &d     = p.second;
		int     start = d.n * d.n * d.n * d.id_local;

		double patch_sum = 0;
		for (int i = 0; i < d.n * d.n; i++) {
			patch_sum += u_view[start + i];
		}

		patch_sum *= d.x_length * d.y_length * d.z_length / (d.n * d.n * d.n);

		sum += patch_sum;
	}
	double retval;
	MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	VecRestoreArray(u, &u_view);
	return retval;
}
double DomainCollection::volume()
{
	double sum = 0;
	for (auto &p : domains) {
		Domain &d = p.second;
		sum += d.x_length * d.y_length * d.z_length;
	}
	double retval;
	MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return retval;
}
PW_explicit<Vec> DomainCollection::getNewDomainVec()
{
	PW<Vec> u;
	VecCreateMPI(MPI_COMM_WORLD, domains.size() * n * n * n, PETSC_DETERMINE, &u);
	return u;
}
