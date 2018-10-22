#ifndef PATCHOPERATOR_H
#define PATCHOPERATOR_H
#include "Domain.h"
#include <petscvec.h>
class PatchOperator
{
	public:
	virtual ~PatchOperator() {}
	virtual void apply(SchurDomain<3> &d, const Vec u, const Vec gamma, Vec f) = 0;
};
#endif
