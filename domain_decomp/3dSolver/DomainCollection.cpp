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
DomainCollection::DomainCollection(OctTree t, int level, int n)
{
	this->n = n;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
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
			shared_ptr<Domain> d_ptr(new Domain());
			Domain &           d = *d_ptr;
			OctNode            n = t.nodes[q.front()];
			q.pop_front();

			d.n        = this->n;
			d.id       = n.id;
			d.x_length = n.x_length;
			d.y_length = n.y_length;
			d.z_length = n.z_length;
			d.x_start  = n.x_start;
			d.y_start  = n.y_start;
			d.z_start  = n.z_start;
			d.child_id = n.child_id;
			if (n.level < level) {
				d.parent_id = n.id;
			} else {
				d.parent_id = n.parent;
			}
			if (d.parent_id != -1) {
				d.oct_on_parent = 0;
				while (t.nodes[n.parent].child_id[d.oct_on_parent] != n.id) {
					d.oct_on_parent++;
				}
			}

			// set and enqueue nbrs
			/*for (Side s : Side::getValues()) {
			    if (n.nbrId(s) != -1) {
			        int id = n.nbrId(s);
			        if (!qed.count(id)) {
			            q.push_back(id);
			            qed.insert(id);
			        }
			        d.getNbrInfoPtr(s) = new NormalNbrInfo(id);
			    }
			}*/
			// set and enqueue nbrs
			for (Side s : Side::getValues()) {
				if (n.nbrId(s) == -1 && n.parent != -1 && t.nodes[n.parent].nbrId(s) != -1) {
					OctNode parent = t.nodes[n.parent];
					OctNode nbr    = t.nodes[parent.nbrId(s)];
					auto    octs   = Octant::getValuesOnSide(s);
					int     quad   = 0;
					while (parent.childId(octs[quad]) != n.id) {
						quad++;
					}
					d.getNbrInfoPtr(s) = new CoarseNbrInfo(nbr.id, quad);
					if (!qed.count(nbr.id)) {
						q.push_back(nbr.id);
						qed.insert(nbr.id);
					}
				} else if (n.level < level && n.nbrId(s) != -1
				           && t.nodes[n.nbrId(s)].hasChildren()) {
					OctNode       nbr  = t.nodes[n.nbrId(s)];
					auto          octs = Octant::getValuesOnSide(s.opposite());
					array<int, 4> nbr_ids;
					for (int i = 0; i < 4; i++) {
						int id     = nbr.childId(octs[i]);
						nbr_ids[i] = id;
						if (!qed.count(id)) {
							q.push_back(id);
							qed.insert(id);
						}
					}
					d.getNbrInfoPtr(s) = new FineNbrInfo(nbr_ids);
				} else if (n.nbrId(s) != -1) {
					int id = n.nbrId(s);
					if (!qed.count(id)) {
						q.push_back(id);
						qed.insert(id);
					}
					d.getNbrInfoPtr(s) = new NormalNbrInfo(id);
				}
			}
			domains[d.id] = d_ptr;
		}
	}
	int num_local_domains = domains.size();
	MPI_Allreduce(&num_local_domains, &num_global_domains, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	reIndex();

	for (auto &p : domains) {
		p.second->setPtrs(domains);
	}
}
DomainCollection::DomainCollection(OctTree t, int n)
{
	this->n = n;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
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
			shared_ptr<Domain> d_ptr(new Domain());
			Domain &           d = *d_ptr;
			OctNode            n = t.nodes[q.front()];
			q.pop_front();

			d.id        = n.id;
			d.n         = this->n;
			d.x_length  = n.x_length;
			d.y_length  = n.y_length;
			d.z_length  = n.z_length;
			d.x_start   = n.x_start;
			d.y_start   = n.y_start;
			d.z_start   = n.z_start;
			d.parent_id = n.parent;
			if (d.parent_id != -1) {
				d.oct_on_parent = 0;
				while (t.nodes[d.parent_id].child_id[d.oct_on_parent] != d.id) {
					d.oct_on_parent++;
				}
			}

			// set and enqueue nbrs
			for (Side s : Side::getValues()) {
				if (n.nbrId(s) == -1 && n.parent != -1 && t.nodes[n.parent].nbrId(s) != -1) {
					OctNode parent = t.nodes[n.parent];
					OctNode nbr    = t.nodes[parent.nbrId(s)];
					auto    octs   = Octant::getValuesOnSide(s);
					int     quad   = 0;
					while (parent.childId(octs[quad]) != n.id) {
						quad++;
					}
					d.getNbrInfoPtr(s) = new CoarseNbrInfo(nbr.id, quad);
					if (!qed.count(nbr.id)) {
						q.push_back(nbr.id);
						qed.insert(nbr.id);
					}
				} else if (n.nbrId(s) != -1 && t.nodes[n.nbrId(s)].hasChildren()) {
					OctNode       nbr  = t.nodes[n.nbrId(s)];
					auto          octs = Octant::getValuesOnSide(s.opposite());
					array<int, 4> nbr_ids;
					for (int i = 0; i < 4; i++) {
						int id     = nbr.childId(octs[i]);
						nbr_ids[i] = id;
						if (!qed.count(id)) {
							q.push_back(id);
							qed.insert(id);
						}
					}
					d.getNbrInfoPtr(s) = new FineNbrInfo(nbr_ids);
				} else if (n.nbrId(s) != -1) {
					int id = n.nbrId(s);
					if (!qed.count(id)) {
						q.push_back(id);
						qed.insert(id);
					}
					d.getNbrInfoPtr(s) = new NormalNbrInfo(id);
				}
			}
			domains[d.id] = d_ptr;
		}
	}
	int num_local_domains = domains.size();
	MPI_Allreduce(&num_local_domains, &num_global_domains, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	reIndex();

	for (auto &p : domains) {
		p.second->setPtrs(domains);
	}
}
DomainCollection::DomainCollection(int d_x, int d_y, int d_z, int n)
{
	this->n    = n;
	auto getID = [&](int x, int y, int z) { return x + y * d_y + z * d_z * d_z; };
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	num_global_domains = d_x * d_y * d_z;
	if (rank == 0) {
		for (int domain_z = 0; domain_z < d_z; domain_z++) {
			for (int domain_y = 0; domain_y < d_y; domain_y++) {
				for (int domain_x = 0; domain_x < d_x; domain_x++) {
					shared_ptr<Domain> d_ptr(new Domain());
					Domain &           ds = *d_ptr;
					ds.n                  = n;
					ds.id                 = getID(domain_x, domain_y, domain_z);
					ds.refine_level       = 1;
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
					domains[ds.id] = d_ptr;
				}
			}
		}
	}
	reIndex();

	for (auto &p : domains) {
		p.second->setPtrs(domains);
	}
}
void DomainCollection::reIndex()
{
	indexDomainsLocal();
}
void DomainCollection::zoltanBalance()
{
	zoltanBalanceDomains();
	for (auto &p : domains) {
		p.second->setPtrs(domains);
	}
	reIndex();
}
void DomainCollection::zoltanBalanceDomains()
{
	struct Zoltan_Struct *zz = Zoltan_Create(MPI_COMM_WORLD);

	// parameters
	Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");  /* Zoltan method: HYPERGRAPH */
	Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); /* Zoltan method: "BLOCK" */
	Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");     /* global ID is 1 integer */
	Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");     /* don't use local IDs */
	Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");      /* we omit object weights */
	Zoltan_Set_Param(zz, "AUTO_MIGRATE", "FALSE");    /* we omit object weights */
	Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");         /* we omit object weights */

	// Query functions
	// Number of Vertices
	auto numObjFn = [](void *data, int *ierr) -> int {
		DomainCollection &dc = *(DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		return dc.domains.size();
	};
	Zoltan_Set_Num_Obj_Fn(zz, numObjFn, this);

	// List of vertices
	auto objListFn
	= [](void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids,
	     ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr) {
		  DomainCollection &dc = *(DomainCollection *) data;
		  *ierr                = ZOLTAN_OK;
		  int pos              = 0;
		  for (auto p : dc.domains) {
			  global_ids[pos] = p.first;
			  pos++;
		  }
	  };
	Zoltan_Set_Obj_List_Fn(zz, objListFn, this);

	// Construct hypergraph
	struct CompressedVertex {
		vector<int> vertices;
		vector<int> ptrs;
		vector<int> edges;
	};
	CompressedVertex graph;
	for (auto &p : domains) {
		Domain &d = *p.second;
		graph.vertices.push_back(d.id);
		graph.ptrs.push_back(graph.edges.size());
		for (Side s : Side::getValues()) {
			int edge_id = -1;
			if (d.hasNbr(s)) {
				switch (d.getNbrType(s)) {
					case NbrType::Normal:
						if (s.isLowerOnAxis()) {
							edge_id = d.getNormalNbrInfo(s).id;
						} else {
							edge_id = d.id;
						}
						break;
					case NbrType::Fine:
						edge_id = d.id;
						break;
					case NbrType::Coarse:
						edge_id = d.getCoarseNbrInfo(s).id;
						break;
				}
			}
			graph.edges.push_back(edge_id * 6 + s.toInt());
		}
	}

	// set graph functions
	Zoltan_Set_HG_Size_CS_Fn(zz,
	                         [](void *data, int *num_lists, int *num_pins, int *format, int *ierr) {
		                         CompressedVertex &graph = *(CompressedVertex *) data;
		                         *ierr                   = ZOLTAN_OK;
		                         *num_lists              = graph.vertices.size();
		                         *num_pins               = graph.edges.size();
		                         *format                 = ZOLTAN_COMPRESSED_VERTEX;
	                         },
	                         &graph);
	Zoltan_Set_HG_CS_Fn(zz,
	                    [](void *data, int num_gid_entries, int num_vtx_edge, int num_pins,
	                       int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr,
	                       ZOLTAN_ID_PTR pin_GID, int *ierr) {
		                    CompressedVertex &graph = *(CompressedVertex *) data;
		                    *ierr                   = ZOLTAN_OK;
		                    copy(graph.vertices.begin(), graph.vertices.end(), vtxedge_GID);
		                    copy(graph.ptrs.begin(), graph.ptrs.end(), vtxedge_ptr);
		                    copy(graph.edges.begin(), graph.edges.end(), pin_GID);
	                    },
	                    &graph);

	Zoltan_Set_Obj_Size_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr) {
		                       DomainCollection &dc = *(DomainCollection *) data;
		                       *ierr                = ZOLTAN_OK;
		                       return dc.domains[*global_id]->serialize(nullptr);
	                       },
	                       this);
	Zoltan_Set_Pack_Obj_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size,
	                          char *buf, int *ierr) {
		                       DomainCollection &dc = *(DomainCollection *) data;
		                       *ierr                = ZOLTAN_OK;
		                       dc.domains[*global_id]->serialize(buf);
		                       dc.domains.erase(*global_id);
	                       },
	                       this);
	Zoltan_Set_Unpack_Obj_Fn(
	zz,
	[](void *data, int num_gid_entries, ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr) {
		DomainCollection &dc = *(DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		dc.domains[*global_id].reset(new Domain());
		dc.domains[*global_id]->deserialize(buf);
	},
	this);
	// Zoltan_Set_Obj_Size_Multi_Fn(zz, DomainZoltanHelper::object_sizes, this);
	// Zoltan_Set_Pack_Obj_Multi_Fn(zz, DomainZoltanHelper::pack_objects, this);
	// Zoltan_Set_Unpack_Obj_Multi_Fn(zz, DomainZoltanHelper::unpack_objects, this);

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

	// update ranks of neighbors before migrating
	for (int i = 0; i < numExport; i++) {
		int curr_id   = exportGlobalIds[i];
		int dest_rank = exportProcs[i];
		domains[curr_id]->updateRank(dest_rank);
	}

	rc = Zoltan_Migrate(zz, numImport, importGlobalIds, importLocalIds, importProcs, importToPart,
	                    numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

	if (rc != ZOLTAN_OK) {
		cerr << "zoltan error\n";
		Zoltan_Destroy(&zz);
		exit(0);
	}
	Zoltan_Destroy(&zz);
#if DD_DEBUG
	cout << "I have " << domains.size() << " domains: ";

	int  prev  = -100;
	bool range = false;
	for (auto &p : domains) {
		int curr = p.second->id;
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
void DomainCollection::zoltanBalanceWithLower(DomainCollection &lower)
{
	struct Zoltan_Struct *zz = Zoltan_Create(MPI_COMM_WORLD);

	// parameters
	Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");  /* Zoltan method: HYPERGRAPH */
	Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); /* Zoltan method: "BLOCK" */
	Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");     /* global ID is 1 integer */
	Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");     /* don't use local IDs */
	Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");      /* we omit object weights */
	Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");     /* we omit object weights */
	Zoltan_Set_Param(zz, "AUTO_MIGRATE", "FALSE");    /* we omit object weights */
	Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");         /* we omit object weights */

	// Query functions
	// Number of Vertices
	struct Levels {
		DomainCollection *upper;
		DomainCollection *lower;
	};
	Levels levels = {this, &lower};
	Zoltan_Set_Num_Obj_Fn(zz,
	                      [](void *data, int *ierr) -> int {
		                      Levels *levels = (Levels *) data;
		                      *ierr          = ZOLTAN_OK;
		                      return levels->upper->domains.size() + levels->lower->domains.size();
	                      },
	                      &levels);

	// List of vertices
	Zoltan_Set_Obj_List_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim,
	                          float *obj_wgts, int *ierr) {
		                       Levels *levels = (Levels *) data;
		                       *ierr          = ZOLTAN_OK;
		                       int pos        = 0;
		                       for (auto p : levels->upper->domains) {
			                       global_ids[pos] = p.first;
			                       obj_wgts[pos]   = 1;
			                       pos++;
		                       }
		                       for (auto p : levels->lower->domains) {
			                       global_ids[pos] = p.first;
			                       obj_wgts[pos]   = 0;
			                       pos++;
		                       }
	                       },
	                       &levels);

	// Construct hypergraph
	struct CompressedVertex {
		vector<int> vertices;
		vector<int> ptrs;
		vector<int> edges;
	};
	CompressedVertex graph;
	// process coarse level
	for (auto &p : domains) {
		Domain &d = *p.second;
		graph.vertices.push_back(d.id);
		graph.ptrs.push_back(graph.edges.size());
		// patch to patch communication
		for (Side s : Side::getValues()) {
			int edge_id = -1;
			if (d.hasNbr(s)) {
				switch (d.getNbrType(s)) {
					case NbrType::Normal:
						if (s.isLowerOnAxis()) {
							edge_id = d.getNormalNbrInfo(s).id;
						} else {
							edge_id = d.id;
						}
						break;
					case NbrType::Fine:
						edge_id = d.id;
						break;
					case NbrType::Coarse:
						edge_id = d.getCoarseNbrInfo(s).id;
						break;
				}
			}
			graph.edges.push_back(edge_id * 6 + s.toInt());
		}
		// level to level communication
		graph.edges.push_back(-d.id - 1);
	}
	// process fine level
	for (auto &p : lower.domains) {
		Domain &d = *p.second;
		graph.vertices.push_back(d.id);
		graph.ptrs.push_back(graph.edges.size());
		graph.edges.push_back(-d.parent_id - 1);
	}

	// set graph functions
	Zoltan_Set_HG_Size_CS_Fn(zz,
	                         [](void *data, int *num_lists, int *num_pins, int *format, int *ierr) {
		                         CompressedVertex &graph = *(CompressedVertex *) data;
		                         *ierr                   = ZOLTAN_OK;
		                         *num_lists              = graph.vertices.size();
		                         *num_pins               = graph.edges.size();
		                         *format                 = ZOLTAN_COMPRESSED_VERTEX;
	                         },
	                         &graph);
	Zoltan_Set_HG_CS_Fn(zz,
	                    [](void *data, int num_gid_entries, int num_vtx_edge, int num_pins,
	                       int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr,
	                       ZOLTAN_ID_PTR pin_GID, int *ierr) {
		                    CompressedVertex &graph = *(CompressedVertex *) data;
		                    *ierr                   = ZOLTAN_OK;
		                    copy(graph.vertices.begin(), graph.vertices.end(), vtxedge_GID);
		                    copy(graph.ptrs.begin(), graph.ptrs.end(), vtxedge_ptr);
		                    copy(graph.edges.begin(), graph.edges.end(), pin_GID);
	                    },
	                    &graph);

	Zoltan_Set_Obj_Size_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr) {
		                       DomainCollection &dc = *(DomainCollection *) data;
		                       *ierr                = ZOLTAN_OK;
		                       return dc.domains.at(*global_id)->serialize(nullptr);
	                       },
	                       this);

	// fixed objects
	Zoltan_Set_Num_Fixed_Obj_Fn(zz,
	                            [](void *data, int *ierr) -> int {
		                            Levels *levels = (Levels *) data;
		                            *ierr          = ZOLTAN_OK;
		                            return levels->lower->domains.size();
	                            },
	                            &levels);

	Zoltan_Set_Fixed_Obj_List_Fn(zz,
	                             [](void *data, int num_fixed_obj, int num_gid_entries,
	                                ZOLTAN_ID_PTR fixed_gids, int *fixed_parts, int *ierr) {
		                             Levels *levels = (Levels *) data;
		                             *ierr          = ZOLTAN_OK;
		                             int rank;
		                             MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		                             int pos = 0;
		                             for (auto p : levels->lower->domains) {
			                             fixed_gids[pos]  = p.first;
			                             fixed_parts[pos] = rank;
			                             pos++;
		                             }
	                             },
	                             &levels);
	// pack and unpack
	// pack and unpack
	Zoltan_Set_Pack_Obj_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size,
	                          char *buf, int *ierr) {
		                       DomainCollection &dc = *(DomainCollection *) data;
		                       *ierr                = ZOLTAN_OK;
		                       dc.domains.at(*global_id)->serialize(buf);
		                       dc.domains.erase(*global_id);
	                       },
	                       this);
	Zoltan_Set_Unpack_Obj_Fn(
	zz,
	[](void *data, int num_gid_entries, ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr) {
		DomainCollection &dc = *(DomainCollection *) data;
		*ierr                = ZOLTAN_OK;
		dc.domains[*global_id].reset(new Domain());
		dc.domains[*global_id]->deserialize(buf);
	},
	this);

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

	// update ranks of neighbors before migrating
	for (int i = 0; i < numExport; i++) {
		int curr_id   = exportGlobalIds[i];
		int dest_rank = exportProcs[i];
		domains.at(curr_id)->updateRank(dest_rank);
	}

	rc = Zoltan_Migrate(zz, numImport, importGlobalIds, importLocalIds, importProcs, importToPart,
	                    numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

	if (rc != ZOLTAN_OK) {
		cerr << "zoltan error\n";
		Zoltan_Destroy(&zz);
		exit(0);
	}
	Zoltan_Destroy(&zz);
#if DD_DEBUG
	cout << "I have " << domains.size() << " domains: ";

	int  prev  = -100;
	bool range = false;
	for (auto &p : domains) {
		int curr = p.second->id;
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

	for (auto &p : domains) {
		p.second->setPtrs(domains);
	}
	reIndex();
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
		p.second->setGlobalNeighborIndexes(rev_map);
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
	int           curr_i = 0;
	vector<int>   map_vec;
	vector<int>   off_proc_map_vec;
	map<int, int> rev_map;
	set<int>      offs;
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
				Domain &d  = *domains[i];
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
		p.second->setLocalNeighborIndexes(rev_map);
	}
	// domain_rev_map          = rev_map;
	domain_map_vec          = map_vec;
	domain_gid_map_vec      = map_vec;
	domain_off_proc_map_vec = off_proc_map_vec;
	indexDomainsGlobal();
}
double DomainCollection::integrate(const Vec u)
{
	double  sum = 0;
	double *u_view;
	VecGetArray(u, &u_view);

	for (auto &p : domains) {
		Domain &d     = *p.second;
		int     start = d.n * d.n * d.n * d.id_local;

		double patch_sum = 0;
		for (int i = 0; i < d.n * d.n * d.n; i++) {
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
		Domain &d = *p.second;
		sum += d.x_length * d.y_length * d.z_length;
	}
	double retval;
	MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return retval;
}
PW_explicit<Vec> DomainCollection::getNewDomainVec() const
{
	PW<Vec> u;
	VecCreateMPI(MPI_COMM_WORLD, domains.size() * n * n * n, PETSC_DETERMINE, &u);
	return u;
}
