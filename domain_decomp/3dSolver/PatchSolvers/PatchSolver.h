#ifndef PATCHSOLVER_H
#define PATCHSOLVER_H
#include "SchurDomain.h"
#include <petscvec.h>
class PatchSolver
{
	public:
	virtual ~PatchSolver() {}
	virtual void addDomain(SchurDomain &d) = 0;
	virtual void domainSolve(std::deque<SchurDomain> &domains, const Vec f, Vec u, const Vec gamma)
	= 0;
};
#endif
