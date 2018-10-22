#include "InterLevelComm.h"
#include <petscao.h>
using namespace std;
using namespace GMG;
InterLevelComm::InterLevelComm(shared_ptr<DomainCollection> coarse_dc,
                               shared_ptr<DomainCollection> fine_dc)
{
	n = coarse_dc->getN();
	set<int> parent_ids;
	for (auto &p : fine_dc->domains) {
		Domain<3> &d = *p.second;
		parent_ids.insert(d.parent_id);
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
		Domain<3> &             d    = *p.second;
		int                     gid  = d.parent_id;
		ILCFineToCoarseMetadata data = {p.second, gid_to_local[gid], gid_to_global[gid]};
		coarse_domains.insert(data);
	}

	PW<IS> dist_is;
	ISCreateBlock(MPI_COMM_SELF, n * n * n, coarse_parent_global_index_map_vec.size(),
	              &coarse_parent_global_index_map_vec[0], PETSC_COPY_VALUES, &dist_is);

	local_vec_size = coarse_parent_global_index_map_vec.size() * n * n * n;

	PW<Vec> u_local = getNewCoarseDistVec();
	PW<Vec> u       = coarse_dc->getNewDomainVec();
	VecScatterCreate(u, dist_is, u_local, nullptr, &scatter);
}
PW_explicit<Vec> InterLevelComm::getNewCoarseDistVec()
{
	PW<Vec> u;
	VecCreateSeq(PETSC_COMM_SELF, local_vec_size, &u);
	return u;
}
PW_explicit<VecScatter> InterLevelComm::getScatter()
{
	return scatter;
}
