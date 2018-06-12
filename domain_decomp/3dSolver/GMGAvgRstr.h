#ifndef GMGAVGRSTR_H
#define GMGAVGRSTR_H
#include "DomainCollection.h"
#include "GMGRestrictor.h"
#include "InterLevelComm.h"
#include <memory>
class GMGAvgRstr : public GMGRestrictor
{
	private:
	std::shared_ptr<DomainCollection> coarse_dc;
	std::shared_ptr<DomainCollection> fine_dc;
	std::shared_ptr<InterLevelComm>   ilc;

	public:
	GMGAvgRstr(std::shared_ptr<DomainCollection> coarse_dc,
	           std::shared_ptr<DomainCollection> fine_dc, std::shared_ptr<InterLevelComm> ilc);
	void restrict(PW<Vec> coarse, PW<Vec> fine);
};
#endif
