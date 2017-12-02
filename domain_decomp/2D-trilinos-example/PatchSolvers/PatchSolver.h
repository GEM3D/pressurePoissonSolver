#ifndef PATCHSOLVER_H
#define PATCHSOLVER_H
#include "Domain.h"
#include "MyTypeDefs.h"
class PatchSolver
{
	public:
	virtual ~PatchSolver() {}
	virtual void addDomain(Domain &d) = 0;
	virtual void solve(Domain &d, const vector_type &f, vector_type &u, const vector_type &gamma)
	= 0;
};
#endif
