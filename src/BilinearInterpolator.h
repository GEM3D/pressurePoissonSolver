#ifndef BILINEARINTERPOLATOR_H
#define BILINEARINTERPOLATOR_H
#include "Interpolator.h"
class BilinearInterpolator : public Interpolator<2>
{
	public:
	void interpolate(SchurDomain<2> &d, const Vec u, Vec interp);
	void interpolate(SchurDomain<2> &d, Side<2> s, int local_index, IfaceType itype, const Vec u,
	                 Vec interp);
};
#endif
