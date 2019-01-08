#ifndef FISHPACKPATCHSOLVER_H
#define FISHPACKPATCHSOLVER_H
#include "PatchSolvers/PatchSolver.h"
class FishpackPatchSolver : public PatchSolver<2>
{
	double lambda = 0;

	public:
	FishpackPatchSolver(double lambda = 0) { this->lambda = lambda; }
	~FishpackPatchSolver() {}
	void addDomain(SchurDomain<2> &d) {}
	void domainSolve(std::deque<SchurDomain<2>> &domains, const Vec f, Vec u, const Vec gamma)
	{
		for (SchurDomain<2> &d : domains) {
			solve(d, f, u, gamma);
		}
	}
	void solve(SchurDomain<2> &d, const Vec f, Vec u, const Vec gamma);
};
#endif
