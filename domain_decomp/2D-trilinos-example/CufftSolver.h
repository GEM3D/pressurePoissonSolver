#ifdef __NVCC__
#ifndef CUFFTSOLVER_H
#define CUFFTSOLVER_H
#include "Domain.h"
#include "Solver.h"
#include <bitset>
#include <cufft.h>
#include <map>
#include <valarray>
struct CDomainKey {
	ulong  neumann = 0;
	double h_x     = 0;
	double h_y     = 0;
	int    n       = 0;

	CDomainKey() {}
	CDomainKey(Domain *d)
	{
		this->neumann = d->ds.neumann.to_ulong();
		this->h_x     = d->h_x;
		this->h_y     = d->h_y;
		this->n       = d->n;
	}
	friend bool operator<(const CDomainKey &l, const CDomainKey &r)
	{
		return std::tie(l.neumann, l.h_x, l.h_y, l.n) < std::tie(r.neumann, r.h_x, r.h_y, r.n);
	}
};

class CufftSolver : public Solver
{
	private:
	static bool compareDomains();
	static std::map<CDomainKey, std::valarray<double>> denoms;
	std::valarray<double> *denom_ptr;
	cufftDoubleReal *(*x_forward)(cufftDoubleReal *, int, int);
	cufftDoubleReal *(*y_forward)(cufftDoubleReal *, int, int);
	cufftDoubleReal *(*x_inverse)(cufftDoubleReal *, int, int);
	cufftDoubleReal *(*y_inverse)(cufftDoubleReal *, int, int);
	static int count;

	public:
	static cufftHandle         n2row;
	static cufftHandle         n2col;
	static cufftHandle         n4row;
	static cufftHandle         n4col;
	static cufftDoubleComplex *w_dst2;
	static cufftDoubleComplex *w_idst2;
	CufftSolver(Domain *d);
	~CufftSolver();
	void solve();
};
#endif
#endif
