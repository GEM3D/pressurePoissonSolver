#ifndef FFTWSOLVER_H
#define FFTWSOLVER_H
#include "Domain.h"
#include "Solver.h"
#include <bitset>
#include <fftw3.h>
#include <map>
#include <valarray>
struct DomainKey {
	ulong  neumann = 0;
	double h_x     = 0;
	double h_y     = 0;
	int    n       = 0;

	DomainKey() {}
	DomainKey(Domain *d)
	{
		this->neumann = d->ds.neumann.to_ulong();
		this->h_x     = d->h_x;
		this->h_y     = d->h_y;
		this->n       = d->n;
	}
	friend bool operator<(const DomainKey &l, const DomainKey &r)
	{
		return std::tie(l.neumann, l.h_x, l.h_y, l.n) < std::tie(r.neumann, r.h_x, r.h_y, r.n);
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
	static std::map<DomainKey, std::valarray<double>> denoms;
	std::valarray<double> *denom_ptr;

	public:
	FftwSolver(Domain *d);
	~FftwSolver();
	void solve();
};
#endif
