#ifndef SOLVER_H
#define SOLVER_H
class Solver;
#include "Domain.h"
class Domain;
class Solver
{
	public:
	virtual ~Solver() {}
	virtual void solve() {}
	protected:
	Domain *d;
};
#endif

