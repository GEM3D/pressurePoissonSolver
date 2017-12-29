#ifndef TRILININTERP_H
#define TRILININTERP_H
#include "Interpolator.h"
class TriLinInterp : public Interpolator
{
	public:
	void interpolate(Domain &d, const Vec u, Vec interp);
	void interpolate(Domain &d, Side s, IfaceType itype, const Vec u, Vec interp);
};
#endif
