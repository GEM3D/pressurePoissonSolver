#ifndef PATCHSOLVER_H
#define PATCHSOLVER_H
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
class PatchSolver
{
	public:
	virtual ~PatchSolver() {}
	void addDomain(DomainSignature &d) {}
	virtual void solve(DomainSignature &d, const vector_type &f, vector_type &u,
	                   const vector_type &gamma)
	{
	}
};
#endif
