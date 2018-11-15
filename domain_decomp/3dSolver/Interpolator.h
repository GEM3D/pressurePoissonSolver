#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include "Iface.h"
#include "SchurDomain.h"
#include "Side.h"
#include <petscvec.h>
template <size_t D>
class Interpolator
{
	public:
	virtual ~Interpolator() {}
	virtual void interpolate(SchurDomain<D> &d, const Vec u, Vec interp) = 0;
	virtual void interpolate(SchurDomain<D> &d, Side<D> s, int local_index, IfaceType itype,
	                         const Vec u, Vec interp)
	= 0;
};
#endif
