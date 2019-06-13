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

#include "P4estDCG.h"

#include <p4est_iterate.h>

using namespace std;

std::shared_ptr<Domain<2>> P4estDCG::getFinestDC()
{
	return domain_list.front();
}
bool P4estDCG::hasCoarserDC()
{
	return curr_level >= 0;
}
/**
 * @brief p4est iterator to get number of levels
 */
static void get_num_levels(p4est_iter_volume_info_t *info, void *user_data)
{
	int &max_level  = *(int *) user_data;
	int  quad_level = info->quad->level;
	if (quad_level > max_level) { max_level = quad_level; }
}
/**
 * @brief p4est iterator to set gids
 */
static void set_ids(p4est_iter_volume_info_t *info, void *user_data)
{
	int *data = (int *) &info->quad->p.user_data;
	data[0]   = *(int *) user_data + info->quadid;
}
/**
 * Constructor
 */
P4estDCG::P4estDCG(p4est_t *p4est, const std::array<int, 2> &ns, IsNeumannFunc<2> inf,
                   BlockMapFunc bmf)
{
	this->ns  = ns;
	this->bmf = bmf;
	this->inf = inf;

	double x_start;
	double y_start;
	bmf(0, 0, 0, x_start, y_start);
	bmf(0, 1, 1, x_scale, y_scale);
	x_scale -= x_start;
	y_scale -= y_start;

	my_p4est = p4est_copy(p4est, false);

	int max_level = 0;
	p4est_iterate_ext(my_p4est, nullptr, &max_level, get_num_levels, nullptr, nullptr, 0);
	int global_max_level = 0;
	MPI_Allreduce(&max_level, &global_max_level, 1, MPI_INT, MPI_MAX, my_p4est->mpicomm);

	int global_id_start;
	MPI_Scan(&my_p4est->local_num_quadrants, &global_id_start, 1, MPI_INT, MPI_SUM,
	         my_p4est->mpicomm);
	global_id_start -= my_p4est->local_num_quadrants;

	// set global ids
	p4est_iterate_ext(my_p4est, nullptr, &global_id_start, set_ids, nullptr, nullptr, 0);

	num_levels = global_max_level + 1;
	curr_level = global_max_level;

	// generate finest DC
	extractLevel();
}
P4estDCG::~P4estDCG()
{
	p4est_destroy(my_p4est);
}
/**
 * @brief iterator to set rank
 */
static void set_ranks(p4est_iter_volume_info_t *info, void *user_data)
{
	int *data = (int *) &info->quad->p.user_data;
	data[1]   = info->p4est->mpirank;
}
struct create_domains_ctx {
	void **                                       ghost_data;
	std::array<int, 2>                            ns;
	std::map<int, std::shared_ptr<PatchInfo<2>>> *dmap;
	P4estDCG::BlockMapFunc                        bmf;
	double                                        x_scale;
	double                                        y_scale;
};
/**
 * @brief iterator to create domain objects
 */
static void create_domains(p4est_iter_volume_info_t *info, void *user_data)
{
	create_domains_ctx &ctx  = *(create_domains_ctx *) user_data;
	auto &              dmap = *ctx.dmap;

	// create domain object
	int                       global_id = info->quad->p.user_int;
	shared_ptr<PatchInfo<2>> &ptr       = dmap[global_id];
	if (ptr == nullptr) { ptr.reset(new PatchInfo<2>); }
	PatchInfo<2> &pinfo = *ptr;

	pinfo.id = global_id;
	pinfo.local_index
	= info->quadid + p4est_tree_array_index(info->p4est->trees, info->treeid)->quadrants_offset;

	// patch dimensions
	pinfo.ns = ctx.ns;

	// cell spacings
	double length     = (double) P4EST_QUADRANT_LEN(info->quad->level) / P4EST_ROOT_LEN;
	pinfo.spacings[0] = length * ctx.x_scale / ctx.ns[0];
	pinfo.spacings[1] = length * ctx.y_scale / ctx.ns[1];

	// coordinates in block
	double x = (double) info->quad->x / P4EST_ROOT_LEN;
	double y = (double) info->quad->y / P4EST_ROOT_LEN;
	ctx.bmf(info->treeid, x, y, pinfo.starts[0], pinfo.starts[1]);

	// set refinement level
	pinfo.refine_level = info->quad->level;
}
/**
 * @brief iterator to set neighbor info
 */
static void link_domains(p4est_iter_face_info_t *info, void *user_data)
{
	create_domains_ctx &ctx  = *(create_domains_ctx *) user_data;
	auto &              dmap = *ctx.dmap;
	if (info->sides.elem_count == 2) {
		p4est_iter_face_side_t side_info1 = ((p4est_iter_face_side_t *) info->sides.array)[0];
		p4est_iter_face_side_t side_info2 = ((p4est_iter_face_side_t *) info->sides.array)[1];

		auto link_side_to_side
		= [&](p4est_iter_face_side_t side_info1, p4est_iter_face_side_t side_info2) {
			  Side<2> side = side_info1.face;
			  if (side_info1.is_hanging) {
				  // coarse nbr
				  int nbr_id, nbr_rank;
				  if (side_info2.is.full.is_ghost) {
					  int *data = (int *) &ctx.ghost_data[side_info2.is.full.quadid];
					  nbr_id    = data[0];
					  nbr_rank  = data[1];
				  } else {
					  nbr_id   = side_info2.is.full.quad->p.user_int;
					  nbr_rank = info->p4est->mpirank;
				  }
				  for (int i = 0; i < 2; i++) {
					  if (!side_info1.is.hanging.is_ghost[i]) {
						  int           id          = side_info1.is.hanging.quad[i]->p.user_int;
						  PatchInfo<2> &pinfo       = *dmap[id];
						  pinfo.getNbrInfoPtr(side) = new CoarseNbrInfo<2>(nbr_id, i);
						  pinfo.getCoarseNbrInfo(side).updateRank(nbr_rank);
					  }
				  }
			  } else if (side_info2.is_hanging) {
				  // fine nbr
				  if (!side_info1.is.full.is_ghost) {
					  int id = side_info1.is.full.quad->p.user_int;

					  std::array<int, 2> nbr_ids;
					  int                ranks[2];
					  for (int i = 0; i < 2; i++) {
						  if (side_info2.is.hanging.is_ghost[i]) {
							  int *data  = (int *) &ctx.ghost_data[side_info2.is.hanging.quadid[i]];
							  nbr_ids[i] = data[0];
							  ranks[i]   = data[1];
						  } else {
							  nbr_ids[i] = side_info2.is.hanging.quad[i]->p.user_int;
							  ranks[i]   = info->p4est->mpirank;
						  }
					  }
					  PatchInfo<2> &pinfo       = *dmap[id];
					  pinfo.getNbrInfoPtr(side) = new FineNbrInfo<2>(nbr_ids);
					  pinfo.getFineNbrInfo(side).updateRank(ranks[0], 0);
					  pinfo.getFineNbrInfo(side).updateRank(ranks[1], 1);
				  }
			  } else {
				  // normal nbr
				  if (!side_info1.is.full.is_ghost) {
					  int id = side_info1.is.full.quad->p.user_int;
					  int nbr_id, nbr_rank;
					  if (side_info2.is.full.is_ghost) {
						  int *data = (int *) &ctx.ghost_data[side_info2.is.full.quadid];
						  nbr_id    = data[0];
						  nbr_rank  = data[1];
					  } else {
						  nbr_id   = side_info2.is.full.quad->p.user_int;
						  nbr_rank = info->p4est->mpirank;
					  }
					  PatchInfo<2> &pinfo       = *dmap[id];
					  pinfo.getNbrInfoPtr(side) = new NormalNbrInfo<2>(nbr_id);
					  pinfo.getNormalNbrInfo(side).updateRank(nbr_rank);
				  }
			  }
		  };

		// 1 to 2
		link_side_to_side(side_info1, side_info2);
		// 2 to 1
		link_side_to_side(side_info2, side_info1);
	}
}

struct coarsen_domains_ctx {
	int                                 level;
	map<int, shared_ptr<PatchInfo<2>>> *dmap;
};
/**
 * @brief iterator to coarsen to next level
 */
static int coarsen(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quads[])
{
	coarsen_domains_ctx &ctx = *(coarsen_domains_ctx *) p4est->user_pointer;
	return quads[0]->level > ctx.level;
}

/**
 * @brief set parent info in finer level
 */
static void coarsen_replace(p4est_t *p4est, p4est_topidx_t which_tree, int num_outgoing,
                            p4est_quadrant_t *outgoing[], int num_incoming,
                            p4est_quadrant_t *incoming[])
{
	coarsen_domains_ctx &ctx  = *(coarsen_domains_ctx *) p4est->user_pointer;
	auto &               dmap = *ctx.dmap;
	incoming[0]->p.user_int   = outgoing[0]->p.user_int;
	for (int i = 0; i < 4; i++) {
		PatchInfo<2> &pinfo  = *dmap[outgoing[i]->p.user_int];
		pinfo.orth_on_parent = i;
		pinfo.parent_id      = incoming[0]->p.user_int;
	}
}
void P4estDCG::extractLevel()
{
	// coarsen previous
	if (curr_level + 1 < num_levels) {
		coarsen_domains_ctx ctx = {curr_level, &domain_list.back()->getPatchInfoMap()};
		my_p4est->user_pointer  = &ctx;
		p4est_coarsen_ext(my_p4est, false, true, coarsen, nullptr, coarsen_replace);
		for (auto p : domain_list.back()->getPatchInfoMap()) {
			PatchInfo<2> &pinfo = *p.second;
			if (!pinfo.hasCoarseParent()) { pinfo.parent_id = pinfo.id; }
		}
		p4est_partition_ext(my_p4est, true, nullptr);
	}
	// set ranks
	p4est_iterate_ext(my_p4est, nullptr, nullptr, set_ranks, nullptr, nullptr, 0);

	// get ids of ghosts
	p4est_ghost_t *ghost = p4est_ghost_new(my_p4est, P4EST_CONNECT_FACE);
	void *         ghost_data[ghost->ghosts.elem_count];
	p4est_ghost_exchange_data(my_p4est, ghost, ghost_data);

	std::map<int, std::shared_ptr<PatchInfo<2>>> new_level;
	create_domains_ctx ctx = {ghost_data, ns, &new_level, bmf, x_scale, y_scale};
	p4est_iterate_ext(my_p4est, ghost, &ctx, create_domains, link_domains, nullptr, 0);

	p4est_ghost_destroy(ghost);

	for (auto p : new_level) {
		p.second->setPtrs(new_level);
	}
	// create DC object
	domain_list.push_back(shared_ptr<Domain<2>>(new Domain<2>(new_level, true)));
	domain_list.back()->setNeumann(inf);

	curr_level--;
}
std::shared_ptr<Domain<2>> P4estDCG::getCoarserDC()
{
	if (curr_level >= 0) {
		extractLevel();
		return domain_list.back();
	} else {
		return nullptr;
	}
}
