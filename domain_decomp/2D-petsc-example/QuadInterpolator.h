#ifndef QUADINTERPOLATOR_H
#define QUADINTERPOLATOR_H
#include "Interpolator.h"
class QuadInterpolator : public Interpolator
{
	public:
	void interpolate(Domain &d, const Vec u, Vec interp);
	void interpolate(Domain &d, Side s, InterpCase icase, const Vec u, Vec interp);
};
#endif
