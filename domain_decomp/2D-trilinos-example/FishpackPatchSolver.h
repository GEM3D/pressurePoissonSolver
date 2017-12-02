#ifndef FISHPACKPATCHSOLVER_H
#define FISHPACKPATCHSOLVER_H
#include "PatchSolver.h"
class FishpackPatchSolver : public PatchSolver
{
	public:
	~FishpackPatchSolver() {}
	void addDomain(Domain &d) {}
	void solve(Domain &d, const vector_type &f, vector_type &u, const vector_type &gamma);
};
#endif
