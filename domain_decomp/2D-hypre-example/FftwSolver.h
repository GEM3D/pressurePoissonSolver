#ifndef FFTWSOLVER_H
#define FFTWSOLVER_H
#include "Solver.h"
#include "Domain.h"
#include <fftw3.h>
class FftwSolver : public Solver
{
	private:
	fftw_plan             plan1;
	fftw_plan             plan2;
	std::valarray<double> f_copy;
	std::valarray<double> tmp;
	std::valarray<double> u;
	std::valarray<double> denom;

	public:
	FftwSolver(Domain *d);
	~FftwSolver();
	void solve();
};
#endif
