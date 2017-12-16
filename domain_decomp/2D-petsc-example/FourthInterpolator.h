#ifndef FOURTHINTERPOLATOR_H
#define FOURTHINTERPOLATOR_H
#include "Interpolator.h"
class FourthInterpolator : public Interpolator
{
	public:
	void interpolate(Domain &d, const Vec u, Vec interp);
	void interpolate(Domain &d, Side s, InterpCase icase, const Vec u, Vec interp);
};
#endif
