/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef THUDNEREGGDCG_H
#define THUDNEREGGDCG_H
#include <Thunderegg/DomainCollectionGenerator.h>
#include <Thunderegg/OctTree.h>
#include <list>
#include <zoltan.h>
template <size_t D> class ThundereggDCG : public DomainCollectionGenerator<D>
{
	private:
	bool neumann;
	/**
	 * @brief finest is stored in front
	 */
	std::list<std::shared_ptr<Domain<D>>> domain;

	int                curr_level;
	int                num_levels;
	Tree<D>            t;
	std::array<int, D> ns;

	void extractLevel();

	using DomainMap = std::map<int, std::shared_ptr<PatchInfo<D>>>;
	void balanceLevel(DomainMap &level);
	void balanceLevelWithLower(DomainMap &level, DomainMap &lower_level);

	public:
	ThundereggDCG(Tree<D> t, std::array<int, D> ns, bool neumann = false);
	~ThundereggDCG() = default;
	std::shared_ptr<Domain<D>> getFinestDC();
	bool                       hasCoarserDC();
	std::shared_ptr<Domain<D>> getCoarserDC();
};
template <size_t D> ThundereggDCG<D>::ThundereggDCG(Tree<D> t, std::array<int, D> ns, bool neumann)
{
	this->t       = t;
	this->ns      = ns;
	this->neumann = neumann;
	num_levels    = t.num_levels;
	curr_level    = num_levels;

	// generate finest DC
	extractLevel();
}
template <size_t D> std::shared_ptr<Domain<D>> ThundereggDCG<D>::getFinestDC()
{
	return domain.front();
}
template <size_t D> std::shared_ptr<Domain<D>> ThundereggDCG<D>::getCoarserDC()
{
	if (curr_level > 0) {
		extractLevel();
		return domain.back();
	} else {
		return nullptr;
	}
}
template <size_t D> bool ThundereggDCG<D>::hasCoarserDC()
{
	return curr_level > 0;
}
template <size_t D> inline void ThundereggDCG<D>::extractLevel()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	DomainMap new_level;
	if (rank == 0) {
		Node<D>         child = *t.levels.at(curr_level);
		std::deque<int> q;
		std::set<int>   qed;
		q.push_back(child.id);
		qed.insert(child.id);

		while (!q.empty()) {
			std::shared_ptr<PatchInfo<D>> d_ptr(new PatchInfo<D>());
			PatchInfo<D> &                pinfo = *d_ptr;
			Node<D>                       n     = t.nodes.at(q.front());
			q.pop_front();

			pinfo.ns = ns;
			pinfo.id = n.id;
			for (size_t i = 0; i < D; i++) {
				pinfo.spacings[i] = n.lengths[i] / ns[i];
			}
			pinfo.starts       = n.starts;
			pinfo.refine_level = n.level;
			if (n.level < curr_level) {
				pinfo.parent_id = n.id;
			} else {
				pinfo.parent_id = n.parent;
				if (pinfo.parent_id != -1) {
					char orth_on_parent = 0;
					while (t.nodes.at(n.parent).child_id[orth_on_parent] != n.id) {
						orth_on_parent++;
					}
					pinfo.orth_on_parent = orth_on_parent;
				}
			}

			// set and enqueue nbrs
			for (Side<D> s : Side<D>::getValues()) {
				if (n.nbrId(s) == -1 && n.parent != -1 && t.nodes.at(n.parent).nbrId(s) != -1) {
					Node<D> parent = t.nodes.at(n.parent);
					Node<D> nbr    = t.nodes.at(parent.nbrId(s));
					auto    octs   = Orthant<D>::getValuesOnSide(s);
					int     quad   = 0;
					while (parent.childId(octs[quad]) != n.id) {
						quad++;
					}
					pinfo.getNbrInfoPtr(s) = new CoarseNbrInfo<D>(nbr.id, quad);
					if (!qed.count(nbr.id)) {
						q.push_back(nbr.id);
						qed.insert(nbr.id);
					}
				} else if (n.level < curr_level && n.nbrId(s) != -1
				           && t.nodes.at(n.nbrId(s)).hasChildren()) {
					Node<D> nbr  = t.nodes.at(n.nbrId(s));
					auto    octs = Orthant<D>::getValuesOnSide(s.opposite());
					std::array<int, Orthant<D>::num_orthants / 2> nbr_ids;
					for (size_t i = 0; i < Orthant<D>::num_orthants / 2; i++) {
						int id     = nbr.childId(octs[i]);
						nbr_ids[i] = id;
						if (!qed.count(id)) {
							q.push_back(id);
							qed.insert(id);
						}
					}
					pinfo.getNbrInfoPtr(s) = new FineNbrInfo<D>(nbr_ids);
				} else if (n.nbrId(s) != -1) {
					int id = n.nbrId(s);
					if (!qed.count(id)) {
						q.push_back(id);
						qed.insert(id);
					}
					pinfo.getNbrInfoPtr(s) = new NormalNbrInfo<D>(id);
				}
			}
			new_level[pinfo.id] = d_ptr;
		}
		for (auto &p : new_level) {
			p.second->setPtrs(new_level);
		}
	}
	// balance
	if (curr_level == num_levels) {
		balanceLevel(new_level);
	} else {
		balanceLevelWithLower(new_level, domain.back()->getPatchInfoMap());
	}
	domain.push_back(std::shared_ptr<Domain<D>>(new Domain<D>(new_level)));
	if (neumann) {
		IsNeumannFunc<D> inf = [](Side<D>, const std::array<double, D> &,
		                          const std::array<double, D> &) { return true; };
		domain.back()->setNeumann(inf);
	}
	curr_level--;
}
template <size_t D> inline void ThundereggDCG<D>::balanceLevel(DomainMap &level)
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
		DomainMap &map = *(DomainMap *) data;
		*ierr          = ZOLTAN_OK;
		return map.size();
	};
	Zoltan_Set_Num_Obj_Fn(zz, numObjFn, &level);

	// List of vertices
	auto objListFn
	= [](void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids,
	     ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr) {
		  DomainMap &map = *(DomainMap *) data;
		  *ierr          = ZOLTAN_OK;
		  int pos        = 0;
		  for (auto p : map) {
			  global_ids[pos] = p.first;
			  pos++;
		  }
	  };
	Zoltan_Set_Obj_List_Fn(zz, objListFn, &level);

	// Construct hypergraph
	struct CompressedVertex {
		std::vector<int> vertices;
		std::vector<int> ptrs;
		std::vector<int> edges;
	};
	CompressedVertex graph;
	for (auto &p : level) {
		PatchInfo<D> &pinfo = *p.second;
		graph.vertices.push_back(pinfo.id);
		graph.ptrs.push_back(graph.edges.size());
		for (Side<D> s : Side<D>::getValues()) {
			int edge_id = -1;
			if (pinfo.hasNbr(s)) {
				switch (pinfo.getNbrType(s)) {
					case NbrType::Normal:
						if (s.isLowerOnAxis()) {
							edge_id = pinfo.getNormalNbrInfo(s).id;
						} else {
							edge_id = pinfo.id;
						}
						break;
					case NbrType::Fine:
						edge_id = pinfo.id;
						break;
					case NbrType::Coarse:
						edge_id = pinfo.getCoarseNbrInfo(s).id;
						break;
				}
			}
			graph.edges.push_back(edge_id * Side<D>::num_sides + s.toInt());
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
		                       DomainMap &map = *(DomainMap *) data;
		                       *ierr          = ZOLTAN_OK;
		                       return map[*global_id]->serialize(nullptr);
	                       },
	                       &level);
	Zoltan_Set_Pack_Obj_Fn(zz,
	                       [](void *data, int num_gid_entries, int num_lid_entries,
	                          ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size,
	                          char *buf, int *ierr) {
		                       DomainMap &map = *(DomainMap *) data;
		                       *ierr          = ZOLTAN_OK;
		                       map[*global_id]->serialize(buf);
		                       map.erase(*global_id);
	                       },
	                       &level);
	Zoltan_Set_Unpack_Obj_Fn(
	zz,
	[](void *data, int num_gid_entries, ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr) {
		DomainMap &map = *(DomainMap *) data;
		*ierr          = ZOLTAN_OK;
		map[*global_id].reset(new PatchInfo<D>());
		map[*global_id]->deserialize(buf);
	},
	&level);
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
		level[curr_id]->updateRank(dest_rank);
	}

	rc = Zoltan_Migrate(zz, numImport, importGlobalIds, importLocalIds, importProcs, importToPart,
	                    numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

	if (rc != ZOLTAN_OK) {
		std::cerr << "zoltan error\n";
		Zoltan_Destroy(&zz);
		exit(0);
	}
	Zoltan_Destroy(&zz);
#if DD_DEBUG
	std::cout << "I have " << level.size() << " domains: ";

	int  prev  = -100;
	bool range = false;
	for (auto &p : level) {
		int curr = p.second->id;
		if (curr != prev + 1 && !range) {
			std::cout << curr << "-";
			range = true;
		} else if (curr != prev + 1 && range) {
			std::cout << prev << " " << curr << "-";
		}
		prev = curr;
	}

	std::cout << prev << "\n";
	std::cout << std::endl;
#endif
	for (auto &p : level) {
		p.second->setPtrs(level);
	}
}
template <size_t D>
inline void ThundereggDCG<D>::balanceLevelWithLower(DomainMap &level, DomainMap &lower_level)
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
		DomainMap *upper;
		DomainMap *lower;
	};
	Levels levels = {&level, &lower_level};
	Zoltan_Set_Num_Obj_Fn(zz,
	                      [](void *data, int *ierr) -> int {
		                      Levels *levels = (Levels *) data;
		                      *ierr          = ZOLTAN_OK;
		                      return levels->upper->size() + levels->lower->size();
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
		                       for (auto p : *levels->upper) {
			                       global_ids[pos] = p.first;
			                       obj_wgts[pos]   = 1;
			                       pos++;
		                       }
		                       for (auto p : *levels->lower) {
			                       global_ids[pos] = p.first;
			                       obj_wgts[pos]   = 0;
			                       pos++;
		                       }
	                       },
	                       &levels);

	// Construct hypergraph
	struct CompressedVertex {
		std::vector<int> vertices;
		std::vector<int> ptrs;
		std::vector<int> edges;
	};
	CompressedVertex graph;
	// process coarse level
	for (auto &p : level) {
		PatchInfo<D> &pinfo = *p.second;
		graph.vertices.push_back(pinfo.id);
		graph.ptrs.push_back(graph.edges.size());
		// patch to patch communication
		for (Side<D> s : Side<D>::getValues()) {
			int edge_id = -1;
			if (pinfo.hasNbr(s)) {
				switch (pinfo.getNbrType(s)) {
					case NbrType::Normal:
						if (s.isLowerOnAxis()) {
							edge_id = pinfo.getNormalNbrInfo(s).id;
						} else {
							edge_id = pinfo.id;
						}
						break;
					case NbrType::Fine:
						edge_id = pinfo.id;
						break;
					case NbrType::Coarse:
						edge_id = pinfo.getCoarseNbrInfo(s).id;
						break;
				}
			}
			graph.edges.push_back(edge_id * Side<D>::num_sides + s.toInt());
		}
		// level to level communication
		graph.edges.push_back(-pinfo.id - 1);
	}
	// process fine level
	for (auto &p : lower_level) {
		PatchInfo<D> &pinfo = *p.second;
		graph.vertices.push_back(pinfo.id);
		graph.ptrs.push_back(graph.edges.size());
		graph.edges.push_back(-pinfo.parent_id - 1);
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
		                       Levels *levels = (Levels *) data;
		                       *ierr          = ZOLTAN_OK;
		                       return levels->upper->at(*global_id)->serialize(nullptr);
	                       },
	                       &levels);

	// fixed objects
	Zoltan_Set_Num_Fixed_Obj_Fn(zz,
	                            [](void *data, int *ierr) -> int {
		                            Levels *levels = (Levels *) data;
		                            *ierr          = ZOLTAN_OK;
		                            return levels->lower->size();
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
		                             for (auto p : *levels->lower) {
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
		                       Levels *levels = (Levels *) data;
		                       *ierr          = ZOLTAN_OK;
		                       levels->upper->at(*global_id)->serialize(buf);
		                       levels->upper->erase(*global_id);
	                       },
	                       &levels);
	Zoltan_Set_Unpack_Obj_Fn(
	zz,
	[](void *data, int num_gid_entries, ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr) {
		Levels *levels = (Levels *) data;
		*ierr          = ZOLTAN_OK;
		(*levels->upper)[*global_id].reset(new PatchInfo<D>());
		(*levels->upper)[*global_id]->deserialize(buf);
	},
	&levels);

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
		levels.upper->at(curr_id)->updateRank(dest_rank);
	}

	rc = Zoltan_Migrate(zz, numImport, importGlobalIds, importLocalIds, importProcs, importToPart,
	                    numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

	if (rc != ZOLTAN_OK) {
		std::cerr << "zoltan error\n";
		Zoltan_Destroy(&zz);
		exit(0);
	}
	Zoltan_Destroy(&zz);
#if DD_DEBUG
	std::cout << "I have " << levels.upper->size() << " domains: ";

	int  prev  = -100;
	bool range = false;
	for (auto &p : *levels.upper) {
		int curr = p.second->id;
		if (curr != prev + 1 && !range) {
			std::cout << curr << "-";
			range = true;
		} else if (curr != prev + 1 && range) {
			std::cout << prev << " " << curr << "-";
		}
		prev = curr;
	}

	std::cout << prev << "\n";
	std::cout << std::endl;
#endif
	for (auto &p : level) {
		p.second->setPtrs(level);
	}
}
extern template class ThundereggDCG<2>;
extern template class ThundereggDCG<3>;
#endif
