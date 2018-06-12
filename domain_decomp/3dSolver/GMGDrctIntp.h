#ifndef GMGDrctIntp_H
#define GMGDrctIntp_H
#include "DomainCollection.h"
#include "GMGInterpolator.h"
#include "InterLevelComm.h"
#include <memory>
class GMGDrctIntp : public GMGInterpolator
{
	private:
	std::shared_ptr<DomainCollection> coarse_dc;
	std::shared_ptr<DomainCollection> fine_dc;
	std::shared_ptr<InterLevelComm>   ilc;

	public:
	GMGDrctIntp(std::shared_ptr<DomainCollection> coarse_dc,
	            std::shared_ptr<DomainCollection> fine_dc, std::shared_ptr<InterLevelComm> ilc);
	void interpolate(PW<Vec> coarse, PW<Vec> fine);
};
#endif
