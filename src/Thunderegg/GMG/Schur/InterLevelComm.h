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

#ifndef InterLevelComm_H
#define InterLevelComm_H
#include "SchurHelper.h"
namespace GMG::Schur
{
/**
 * @brief Structure that wraps some extra meta-data around a Domain object.
 */
template <size_t D> struct ILCFineToCoarseMetadata {
	/**
	 * @brief The domain that this meta-data corresponds to.
	 */
	std::shared_ptr<IfaceSet<D>> ifs;
	/**
	 * @brief the local block-index of the parent domain in the scattered coarse vector.
	 */
	int local_index;
	/**
	 * @brief the global block-index of the parent domain in the scattered coarse vector.
	 */
	int global_index;
	/**
	 * @brief less than operator so that this struct can be placed in set container.
	 */
	bool operator<(const ILCFineToCoarseMetadata &other) const
	{
		return *d < *other.d;
	}
};
/**
 * @brief Creates a mapping from fine to coarse levels.
 */
template <size_t D> class InterLevelComm
{
	private:
	/**
	 * @brief The number of cells in each direction on a patch
	 */
	int n;
	/**
	 * @brief the number of elements in the local distributed vector.
	 */
	int local_vec_size;
	/**
	 * @brief a set of domains on the finer level that contains extra meta-data on the mapping from
	 * finer level vectors to coarser level vectors.
	 */
	std::set<ILCFineToCoarseMetadata<D>> coarse_domains;
	/**
	 * @brief The PETSc VecScatter object that scatters the values of the coarse vector to each
	 * process that needs them for interpolation / restriction.
	 */
	PW<VecScatter> scatter;

	public:
	/**
	 * @brief Create a new InterLevelComm object.
	 *
	 * @param coarse_dc the coarser DomainCollection.
	 * @param fine_dc the finer DomainCollection.
	 */
	InterLevelComm(std::shared_ptr<SchurHelper<D>> coarse_sh,
	               std::shared_ptr<SchurHelper<D>> fine_sh);

	/**
	 * @brief Allocate a new vector for which coarse values values will be scattered into.
	 *
	 * @return the newly allocated vector.
	 */
	PW_explicit<Vec> getNewCoarseDistVec();

	/**
	 * @brief Get a PETSc VecScatter object that will scatter values from the coarse vector to the
	 * processes that need them for interpolation / restriction.
	 *
	 * @return the VecScatter object
	 */
	PW_explicit<VecScatter> getScatter();

	/**
	 * @brief get a set of domains on the finer level that contains meta-data for how the values in
	 * the fine vector map to the values in the coarse vector.
	 *
	 * @return set of fine domains
	 */
	std::set<ILCFineToCoarseMetadata<D>> getFineDomains()
	{
		return coarse_domains;
	}
};
template <size_t D>
inline InterLevelComm<D>::InterLevelComm(std::shared_ptr<SchurHelper<D>> coarse_sh,
                                         std::shared_ptr<SchurHelper<D>> fine_sh)
{
	using namespace std;
	n = coarse_sh->getN();
	set<int> parent_ids;
	for (auto &p : fine_sh->getIfaces()) {
		IfaceSet<D> &ifs = *p.second;
        if(ifs.parent_id!=-1){
		parent_ids.insert(d.parent_id);
        }
	}
	vector<int> coarse_parent_global_index_map_vec(parent_ids.begin(), parent_ids.end());
	vector<int> coarse_parent_gid_map_vec = coarse_parent_global_index_map_vec;
	// get global indexes for parent domains
	PW<AO> ao;
	AOCreateMapping(MPI_COMM_WORLD, coarse_dc->domain_gid_map_vec.size(),
	                &coarse_dc->domain_gid_map_vec[0], &coarse_dc->domain_map_vec[0], &ao);
	AOApplicationToPetsc(ao, coarse_parent_global_index_map_vec.size(),
	                     &coarse_parent_global_index_map_vec[0]);

	// set index info
	map<int, int> gid_to_local;
	map<int, int> gid_to_global;

	for (size_t i = 0; i < coarse_parent_gid_map_vec.size(); i++) {
		int gid            = coarse_parent_gid_map_vec[i];
		gid_to_local[gid]  = i;
		gid_to_global[gid] = coarse_parent_global_index_map_vec[i];
	}
	for (auto &p : fine_dc->domains) {
		Domain<D> &                d    = *p.second;
		int                        gid  = d.parent_id;
		ILCFineToCoarseMetadata<D> data = {p.second, gid_to_local[gid], gid_to_global[gid]};
		coarse_domains.insert(data);
	}

	PW<IS> dist_is;
	ISCreateBlock(MPI_COMM_SELF, (int) pow(n, D-1), coarse_parent_global_index_map_vec.size(),
	              &coarse_parent_global_index_map_vec[0], PETSC_COPY_VALUES, &dist_is);

	local_vec_size = coarse_parent_global_index_map_vec.size() * pow(n, D-1);

	PW<Vec> u_local = getNewCoarseDistVec();
	PW<Vec> u       = coarse_dc->getNewDomainVec();
	VecScatterCreate(u, dist_is, u_local, nullptr, &scatter);
}

template <size_t D> inline PW_explicit<Vec> InterLevelComm<D>::getNewCoarseDistVec()
{
	PW<Vec> u;
	VecCreateSeq(PETSC_COMM_SELF, local_vec_size, &u);
	return u;
}
template <size_t D> inline PW_explicit<VecScatter> InterLevelComm<D>::getScatter()
{
	return scatter;
}
} // namespace GMG
#endif
