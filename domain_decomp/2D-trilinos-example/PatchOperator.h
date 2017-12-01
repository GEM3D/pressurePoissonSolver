#ifndef PATCHOPERATOR_H
#define PATCHOPERATOR_H
#include "DomainSignatureCollection.h"
class PatchOperator
{
	public:
	virtual ~PatchOperator() {}
	virtual void apply(DomainSignature &d, const vector_type &u, const vector_type &gamma,
	                   vector_type &f)
	= 0;
};
#endif
