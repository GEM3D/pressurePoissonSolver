#ifndef FISHPACKSOLVER_H
#define FISHPACKSOLVER_H
#include "Solver.h"
#include "Domain.h"
#include <fftw3.h>
class FishpackSolver : public Solver
{
	public:
	FishpackSolver(Domain *d);
	~FishpackSolver() {}
	void solve();
};
#endif
