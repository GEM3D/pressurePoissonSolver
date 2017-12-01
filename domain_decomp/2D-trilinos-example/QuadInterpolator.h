#ifndef QUADINTERPOLATOR_H
#define QUADINTERPOLATOR_H
#include "Interpolator.h"
class QuadInterpolator : public Interpolator
{
	private:
	int n;

	public:
	QuadInterpolator(int n);
	void interpolate(DomainSignature &d, const vector_type &u, vector_type &interp);
	void interpolate(DomainSignature &d, Side s, InterpCase icase, const vector_type &u,
	                 vector_type &interp);
};
#endif
