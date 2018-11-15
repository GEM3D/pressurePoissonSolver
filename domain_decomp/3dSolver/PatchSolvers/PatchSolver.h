#ifndef PATCHSOLVER_H
#define PATCHSOLVER_H
#include "SchurDomain.h"
#include <petscvec.h>
template <size_t D>
class PatchSolver
{
	public:
	virtual ~PatchSolver() {}
	virtual void addDomain(SchurDomain<D> &d) = 0;
	virtual void domainSolve(std::deque<SchurDomain<D>> &domains, const Vec f, Vec u, const Vec gamma)
	= 0;
};
#endif
