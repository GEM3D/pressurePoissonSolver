#ifndef QUADINTERPOLATOR_H
#define QUADINTERPOLATOR_H
#include "Interpolator.h"
class QuadInterpolator : public Interpolator<2>
{
	public:
	void interpolate(Domain<2> &d, const Vec u, Vec interp);
	void interpolate(Domain<2> &d, Side<2> s, InterpCase icase, const Vec u, Vec interp);
};
#endif
