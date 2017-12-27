#ifndef PATCHSOLVER_H
#define PATCHSOLVER_H
#include "Domain.h"
#include <petscvec.h>
class PatchSolver
{
	public:
	virtual ~PatchSolver() {}
	virtual void addDomain(Domain &d) = 0;
	virtual void solve(Domain &d, const Vec f, Vec u, const Vec gamma) = 0;
};
#endif
