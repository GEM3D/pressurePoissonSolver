#ifndef FFTWSOLVER_H
#define FFTWSOLVER_H
#include "Domain.h"
#include "Solver.h"
#include <bitset>
#include <fftw3.h>
#include <map>
#include <valarray>
class DomainCmp
{
	public:
	bool operator()(Domain *const &lhs, Domain *const &rhs) const
	{
		ulong l = lhs->neumann_sides.to_ulong();
		ulong r = rhs->neumann_sides.to_ulong();
		return std::tie(l, lhs->h_x, lhs->h_y, lhs->n) < std::tie(r, rhs->h_x, rhs->h_y, rhs->n);
	}
};
class FftwSolver : public Solver
{
	private:
	static bool                  compareDomains();
	fftw_plan                    plan1;
	fftw_plan                    plan2;
	static bool                  initialized;
	static std::valarray<double> f_copy;
	static std::valarray<double> tmp;
	static std::valarray<double> u;
	static std::map<Domain *, std::valarray<double>, DomainCmp> denoms;
	std::valarray<double> *denom_ptr;

	public:
	FftwSolver(Domain *d);
	~FftwSolver();
	void solve();
};
#endif
