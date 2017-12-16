#ifndef QUADINTERPOLATOR_H
#define QUADINTERPOLATOR_H
#include "Interpolator.h"
class QuadInterpolator : public Interpolator
{
	public:
	void interpolate(Domain &d, const vector_type &u, vector_type &interp);
	void interpolate(Domain &d, Side s, InterpCase icase, const vector_type &u,
	                 vector_type &interp);
};
#endif
