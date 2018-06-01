#ifndef GMGAVGRSTR_H
#define GMGAVGRSTR_H
#include "DomainCollection.h"
#include "GMGRestrictor.h"
#include <memory>
class GMGAvgRstr : public GMGRestrictor
{
	private:
	DomainCollection coarse_dc;
	DomainCollection fine_dc;

	public:
	GMGAvgRstr(DomainCollection &coarse_dc, DomainCollection &fine_dc);
	void restrict(PW<Vec> coarse, PW<Vec> fine);
};
#endif
