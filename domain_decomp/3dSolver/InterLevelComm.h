#ifndef InterLevelComm_H
#define InterLevelComm_H
#include "DomainCollection.h"
struct ILCFineToCoarseMetadata {
	std::shared_ptr<Domain> d;
	int                     local_index;
	int                     global_index;
	bool operator<(const ILCFineToCoarseMetadata &other) const { return *d < *other.d; }
};
class InterLevelComm
{
	private:
	int                               n;
	int                               local_vec_size;
	std::set<ILCFineToCoarseMetadata> coarse_domains;
	PW<VecScatter>                    scatter;

	public:
	InterLevelComm(std::shared_ptr<DomainCollection> coarse_dc,
	               std::shared_ptr<DomainCollection> fine_dc);

	PW_explicit<Vec>                  getNewCoarseDistVec();
	PW_explicit<VecScatter>           getScatter();
	std::set<ILCFineToCoarseMetadata> getFineDomains() { return coarse_domains; }
};
#endif
