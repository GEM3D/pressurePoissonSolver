#ifndef FISHPACKPATCHSOLVER_H
#define FISHPACKPATCHSOLVER_H
#include "PatchSolvers/PatchSolver.h"
class FishpackPatchSolver : public PatchSolver
{
	double lambda = 0;

	public:
	FishpackPatchSolver(double lambda = 0) { this->lambda = lambda; }
	~FishpackPatchSolver() {}
	void addDomain(Domain &d) {}
	void solve(Domain &d, const Vec f, Vec u, const Vec gamma);
};
#endif
