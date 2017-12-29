#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include "Domain.h"
#include "Iface.h"
#include "Side.h"
#include <petscvec.h>
class Interpolator
{
	public:
	virtual ~Interpolator() {}
	virtual void interpolate(Domain &d, const Vec u, Vec interp) = 0;
	virtual void interpolate(Domain &d, Side s, IfaceType itype, const Vec u, Vec interp) = 0;
};
#endif
