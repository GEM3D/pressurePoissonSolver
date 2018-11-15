#ifndef TRILININTERP_H
#define TRILININTERP_H
#include "Interpolator.h"
class TriLinInterp : public Interpolator<3>
{
	public:
	void interpolate(SchurDomain<3> &d, const Vec u, Vec interp);
	void interpolate(SchurDomain<3> &d, Side<3> s, int local_index, IfaceType itype, const Vec u,
	                 Vec interp);
};
#endif
