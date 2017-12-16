#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include "Domain.h"
#include "InterpCase.h"
#include "MyTypeDefs.h"
#include "Side.h"
class Interpolator
{
	public:
	virtual ~Interpolator() {}
	virtual void interpolate(Domain &d, const vector_type &u, vector_type &interp) = 0;
	virtual void interpolate(Domain &d, Side s, InterpCase icase, const vector_type &u,
	                         vector_type &interp)
	= 0;
};
#endif
