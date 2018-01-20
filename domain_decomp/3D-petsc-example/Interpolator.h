#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include "Iface.h"
#include "SchurDomain.h"
#include "Side.h"
#include <petscvec.h>
class Interpolator
{
	public:
	virtual ~Interpolator() {}
	virtual void interpolate(SchurDomain &d, const Vec u, Vec interp) = 0;
	virtual void interpolate(SchurDomain &d, Side s, int local_index, IfaceType itype, const Vec u,
	                         Vec interp)
	= 0;
};
#endif
