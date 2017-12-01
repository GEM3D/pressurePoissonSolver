#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include "DomainSignatureCollection.h"
#include "InterpCase.h"
#include "Side.h"
class Interpolator
{
	public:
	virtual ~Interpolator() {}
	virtual void interpolate(DomainSignature &d, const vector_type &u, vector_type &interp) = 0;
	virtual void interpolate(DomainSignature &d, Side s, InterpCase icase, const vector_type &u,
	                         vector_type &interp)
	= 0;
};
#endif
