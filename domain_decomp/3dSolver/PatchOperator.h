#ifndef PATCHOPERATOR_H
#define PATCHOPERATOR_H
#include "SchurDomain.h"
#include <petscvec.h>
template <size_t D> class PatchOperator
{
	public:
	virtual ~PatchOperator() {}
	virtual void apply(SchurDomain<D> &d, const Vec u, const Vec gamma, Vec f) = 0;
};
#endif
