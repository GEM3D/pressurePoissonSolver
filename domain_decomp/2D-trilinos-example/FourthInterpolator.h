#ifndef FOURTHINTERPOLATOR_H
#define FOURTHINTERPOLATOR_H
#include "Interpolator.h"
class FourthInterpolator : public Interpolator
{
	public:
	void interpolate(Domain &d, const vector_type &u, vector_type &interp);
	void interpolate(Domain &d, Side s, InterpCase icase, const vector_type &u,
	                 vector_type &interp);
};
#endif
