#ifndef GMGDrctIntp_H
#define GMGDrctIntp_H
#include "DomainCollection.h"
#include "GMGInterpolator.h"
#include <memory>
class GMGDrctIntp : public GMGInterpolator
{
	private:
	DomainCollection coarse_dc;
	DomainCollection fine_dc;

	public:
	GMGDrctIntp(DomainCollection &coarse_dc, DomainCollection &fine_dc);
	void interpolate(PW<Vec> coarse, PW<Vec> fine);
};
#endif
