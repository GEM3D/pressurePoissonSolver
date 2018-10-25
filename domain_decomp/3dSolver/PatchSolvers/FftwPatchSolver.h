#ifndef FFTWPATCHSOLVER_H
#define FFTWPATCHSOLVER_H
#include "DomainCollection.h"
#include "PatchSolvers/PatchSolver.h"
#include <bitset>
#include <fftw3.h>
#include <map>
#include <valarray>
#ifndef DOMAINK
#define DOMAINK
struct DomainK {
	ulong  neumann = 0;
	double h_x     = 0;
	double h_y     = 0;

	DomainK() {}
	DomainK(const SchurDomain<3> &d)
	{
		this->neumann = d.neumann.to_ulong();
		this->h_x     = d.domain.lengths[0];
		this->h_y     = d.domain.lengths[1];
	}
	friend bool operator<(const DomainK &l, const DomainK &r)
	{
		return std::tie(l.neumann, l.h_x, l.h_y) < std::tie(r.neumann, r.h_x, r.h_y);
	}
};
#endif

class FftwPatchSolver : public PatchSolver
{
	private:
	int                                      n;
	bool                                     initialized = false;
	static bool                              compareDomains();
	double                                   lambda;
	std::map<DomainK, fftw_plan>             plan1;
	std::map<DomainK, fftw_plan>             plan2;
	std::valarray<double>                    f_copy;
	std::valarray<double>                    tmp;
	std::valarray<double>                    sol;
	std::map<DomainK, std::valarray<double>> denoms;

	public:
	FftwPatchSolver(DomainCollection &dsc, double lambda = 0);
	~FftwPatchSolver();
	void solve(SchurDomain<3> &d, const Vec f, Vec u, const Vec gamma);
	void domainSolve(std::deque<SchurDomain<3>> &domains, const Vec f, Vec u, const Vec gamma)
	{
		for (SchurDomain<3> &d : domains) {
			solve(d, f, u, gamma);
		}
	}
	void addDomain(SchurDomain<3> &d);
};
#endif
