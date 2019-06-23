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

#ifndef THUNDEREGG_DOMAIN_H
#define THUNDEREGG_DOMAIN_H
#include <Thunderegg/InterpCase.h>
#include <Thunderegg/PW.h>
#include <Thunderegg/PatchInfo.h>
#include <Thunderegg/PetscVector.h>
#include <cmath>
#include <deque>
#include <map>
#include <memory>
#include <petscao.h>
#include <petscvec.h>
#include <set>
#include <string>
#include <vector>
/**
 * @brief Represents the domain of the problem.
 *
 * This class mainly manages a set of patches that makes up the domain. It is responsible for
 * setting up the indexing of the domains, which is used in the rest of the Thunderegg library.
 *
 * @tparam D the number of Cartesian dimensions
 */
template <size_t D> class Domain
{
	private:
	/**
	 * @brief The number of cells in each direction
	 */
	std::array<int, D> ns;
	/**
	 * @brief The number of cells in a patch
	 */
	int num_cells_in_patch;
	/**
	 * @brief Map that goes form patch's id to the PatchInfo pointer
	 */
	std::map<int, std::shared_ptr<PatchInfo<D>>> pinfo_id_map;
	/**
	 * @brief Vector of PatchInfo pointers where index in the vector corresponds to the patch's
	 * local index
	 */
	std::vector<std::shared_ptr<PatchInfo<D>>> pinfo_vector;
	/**
	 * @brief The global number of patches
	 */
	int global_num_patches = 1;
	/**
	 * @brief Vector of patch ids. Index corresponds to the patch's local index.
	 */
	std::vector<int> patch_id_map_vec;
	/**
	 * @brief Vector of patch global_indexes. Index corresponds to the patch's global index.
	 */
	std::vector<int> patch_global_index_map_vec;
	/**
	 * @brief Vector of ghost patch global_indexes. Index corresponds to the patch's local index in
	 * the ghost vector.
	 */
	std::vector<int> patch_global_index_map_vec_off_proc;
	/**
	 * @brief mpi rank
	 */
	int rank;

	/**
	 * @brief Give the patches local indexes.
	 *
	 * @param local_id_set true if local indexes are set by user.
	 */
	void indexDomainsLocal(bool local_id_set = false);
	/**
	 * @brief Give the patches global indexes
	 *
	 * @param global_id_set true if global indexes are set by user.
	 */
	void indexDomainsGlobal(bool global_id_set = false);

	/**
	 * @brief Give patches new global and local indexes
	 *
	 * @param local_id_set true if local indexes are set by user
	 * @param global_id_set true if global indexes are set by user
	 */
	void reIndex(bool local_id_set, bool global_id_set)
	{
		indexDomainsLocal(local_id_set);
		if (!global_id_set) { indexDomainsGlobal(); }
	}

	public:
	/**
	 * @brief Construct a new Domain object with a given PatchInfo map.
	 *
	 * @param pinfo_map map that goes from patch id to the PatchInfo pointer
	 * @param local_id_set true if local indexes are set by user
	 * @param global_id_set true if global indexes are set by user
	 */
	Domain(std::map<int, std::shared_ptr<PatchInfo<D>>> pinfo_map, bool local_id_set = false,
	       bool global_id_set = false)
	{
		this->ns = pinfo_map.begin()->second->ns;

		num_cells_in_patch = 1;
		for (size_t i = 0; i < D; i++) {
			num_cells_in_patch *= ns[i];
		}

		pinfo_id_map = pinfo_map;

		int num_local_domains = pinfo_id_map.size();
		MPI_Allreduce(&num_local_domains, &global_num_patches, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		reIndex(local_id_set, global_id_set);

		for (auto &p : pinfo_id_map) {
			p.second->setPtrs(pinfo_id_map);
		}
	}
	/**
	 * @brief Get a vector of PatchInfo pointers where index in the vector corresponds to the
	 * patch's local index
	 */
	std::vector<std::shared_ptr<PatchInfo<D>>> &getPatchInfoVector()
	{
		return pinfo_vector;
	}
	/**
	 * @brief Get map that goes form patch's id to the PatchInfo pointer
	 */
	std::map<int, std::shared_ptr<PatchInfo<D>>> &getPatchInfoMap()
	{
		return pinfo_id_map;
	}
	/**
	 * @brief Get a vector of patch ids. Index in vector corresponds to the patch's local index.
	 */
	const std::vector<int> &getIdMapVec() const
	{
		return patch_id_map_vec;
	}
	/**
	 * @brief Get a vector of patch global indexes. Index in vector corresponds to the patch's
	 * global index.
	 */
	const std::vector<int> &getGlobalIndexMapVec() const
	{
		return patch_global_index_map_vec;
	}
	/**
	 * @brief Get the number of cells in each direction
	 *
	 */
	const std::array<int, D> &getNs() const
	{
		return ns;
	}
	/**
	 * @brief Set the neumann boundary conditions
	 *
	 * @param inf the function that determines boundary conditions
	 */
	void setNeumann(IsNeumannFunc<D> inf)
	{
		for (auto &p : pinfo_id_map) {
			p.second->setNeumann(inf);
		}
	}
	std::shared_ptr<PetscVector<D>> getNewDomainVec() const
	{
		Vec u;
		VecCreateMPI(MPI_COMM_WORLD, pinfo_id_map.size() * num_cells_in_patch, PETSC_DETERMINE, &u);
		return std::shared_ptr<PetscVector<D>>(new PetscVector<D>(u, ns));
	}
	/**
	 * @brief Get the number of global patches
	 */
	int getNumGlobalPatches() const
	{
		return global_num_patches;
	}
	/**
	 * @brief Get the number of local patches
	 */
	int getNumLocalPatches() const
	{
		return pinfo_vector.size();
	}
	/**
	 * @brief get the number of global cells
	 */
	int getNumGlobalCells() const
	{
		return global_num_patches * num_cells_in_patch;
	}
	/**
	 * @brief Get get the number of local cells
	 */
	int getNumLocalCells() const
	{
		return pinfo_id_map.size() * num_cells_in_patch;
	}
	/**
	 * @brief Get the number of cells in a patch
	 */
	int getNumCellsInPatch() const
	{
		return num_cells_in_patch;
	}
	/**
	 * @brief Get the volume of the domain.
	 *
	 * For 2D, this will be the area.
	 */
	double volume() const
	{
		double sum = 0;
		for (auto &p : pinfo_id_map) {
			PatchInfo<D> &d         = *p.second;
			double        patch_vol = 1;
			for (size_t i = 0; i < D; i++) {
				patch_vol *= d.spacings[i] * d.ns[i];
			}
			sum += patch_vol;
		}
		double retval;
		MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return retval;
	}
	/**
	 * @brief Integrate a vector over the domain.
	 *
	 * @param u the vector
	 * @return double the result of the integral
	 */
	double integrate(std::shared_ptr<const Vector<D>> u) const
	{
		double sum = 0;

		for (auto &p : pinfo_id_map) {
			PatchInfo<D> &     d      = *p.second;
			const LocalData<D> u_data = u->getLocalData(d.local_index);

			double patch_sum = 0;
			nested_loop<D>(u_data.getStart(), u_data.getEnd(),
			               [&](std::array<int, D> coord) { patch_sum += u_data[coord]; });

			for (size_t i = 0; i < D; i++) {
				patch_sum *= d.spacings[i];
			}
			sum += patch_sum;
		}
		double retval;
		MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return retval;
	}
};

template <size_t D> void Domain<D>::indexDomainsLocal(bool local_id_set)
{
	std::vector<int>   map_vec;
	std::vector<int>   off_proc_map_vec;
	std::map<int, int> rev_map;
	int                curr_i = local_id_set ? pinfo_id_map.size() : 0;
	if (local_id_set) {
		std::set<int> offs;
		if (!pinfo_id_map.empty()) {
			std::set<int> todo;
			map_vec.resize(pinfo_id_map.size());
			for (auto &p : pinfo_id_map) {
				todo.insert(p.first);
				PatchInfo<D> &d        = *p.second;
				map_vec[d.local_index] = d.id;
			}
			std::set<int> enqueued;
			while (!todo.empty()) {
				std::deque<int> queue;
				queue.push_back(*todo.begin());
				enqueued.insert(*todo.begin());
				while (!queue.empty()) {
					int i = queue.front();
					todo.erase(i);
					queue.pop_front();
					PatchInfo<D> &d = *pinfo_id_map[i];
					rev_map[i]      = d.local_index;
					for (int i : d.getNbrIds()) {
						if (!enqueued.count(i)) {
							enqueued.insert(i);
							if (pinfo_id_map.count(i)) {
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
	} else {
		std::set<int> offs;
		if (!pinfo_id_map.empty()) {
			std::set<int> todo;
			for (auto &p : pinfo_id_map) {
				todo.insert(p.first);
			}
			std::set<int> enqueued;
			while (!todo.empty()) {
				std::deque<int> queue;
				queue.push_back(*todo.begin());
				enqueued.insert(*todo.begin());
				while (!queue.empty()) {
					int i = queue.front();
					todo.erase(i);
					queue.pop_front();
					map_vec.push_back(i);
					PatchInfo<D> &d = *pinfo_id_map[i];
					rev_map[i]      = curr_i;
					d.local_index   = curr_i;
					curr_i++;
					for (int i : d.getNbrIds()) {
						if (!enqueued.count(i)) {
							enqueued.insert(i);
							if (pinfo_id_map.count(i)) {
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
	}
	// map off proc
	for (int i : off_proc_map_vec) {
		rev_map[i] = curr_i;
		curr_i++;
	}
	pinfo_vector.resize(pinfo_id_map.size());
	for (auto &p : pinfo_id_map) {
		p.second->setLocalNeighborIndexes(rev_map);
		pinfo_vector[p.second->local_index] = p.second;
	}
	// domain_rev_map          = rev_map;
	patch_global_index_map_vec          = map_vec;
	patch_id_map_vec                    = map_vec;
	patch_global_index_map_vec_off_proc = off_proc_map_vec;
}
template <size_t D> void Domain<D>::indexDomainsGlobal(bool global_id_set)
{
	// global indices are going to be sequentially increasing with rank
	int local_size = pinfo_id_map.size();
	int start_i;
	MPI_Scan(&local_size, &start_i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	start_i -= local_size;
	std::vector<int> new_global(local_size);
	iota(new_global.begin(), new_global.end(), start_i);

	// create map for gids
	PW<AO> ao;
	AOCreateMapping(MPI_COMM_WORLD, local_size, &patch_global_index_map_vec[0], &new_global[0],
	                &ao);

	// get global indices that we want to recieve for dest vector
	std::vector<int> inds = patch_global_index_map_vec;
	for (int i : patch_global_index_map_vec_off_proc) {
		inds.push_back(i);
	}

	// get new global indices
	AOApplicationToPetsc(ao, inds.size(), &inds[0]);
	std::map<int, int> rev_map;
	for (size_t i = 0; i < inds.size(); i++) {
		rev_map[i] = inds[i];
	}

	for (auto &p : pinfo_id_map) {
		p.second->setGlobalNeighborIndexes(rev_map);
	}
	for (size_t i = 0; i < patch_global_index_map_vec.size(); i++) {
		patch_global_index_map_vec[i] = inds[i];
	}
	for (size_t i = 0; i < patch_global_index_map_vec_off_proc.size(); i++) {
		patch_global_index_map_vec_off_proc[i] = inds[patch_global_index_map_vec.size() + i];
	}
}
template <size_t D> class DomainVG : public VectorGenerator<D>
{
	private:
	std::shared_ptr<Domain<D>> dc;

	public:
	DomainVG(std::shared_ptr<Domain<D>> dc)
	{
		this->dc = dc;
	}
	std::shared_ptr<Vector<D>> getNewVector()
	{
		return dc->getNewDomainVec();
	}
};
extern template class Domain<2>;
extern template class Domain<3>;
#endif
