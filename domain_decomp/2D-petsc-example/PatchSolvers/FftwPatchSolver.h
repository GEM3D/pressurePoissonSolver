#ifndef FFTWPATCHSOLVER_H
#define FFTWPATCHSOLVER_H
#include "DomainCollection.h"
#include "PatchSolvers/PatchSolver.h"
#include <bitset>
#include <fftw3.h>
#include <map>
#include <valarray>
struct DomainK {
	ulong  neumann = 0;
	double h_x     = 0;
	double h_y     = 0;

	DomainK() {}
	DomainK(const Domain &d)
	{
		this->neumann = d.neumann.to_ulong();
		this->h_x     = d.x_length;
		this->h_y     = d.y_length;
	}
	friend bool operator<(const DomainK &l, const DomainK &r)
	{
		return std::tie(l.neumann, l.h_x, l.h_y) < std::tie(r.neumann, r.h_x, r.h_y);
	}
};

class FftwPatchSolver : public PatchSolver
{
	private:
	int         n;
	bool        initialized = false;
	static bool compareDomains();
    double lambda;
	std::map<DomainK, fftw_plan> plan1;
	std::map<DomainK, fftw_plan> plan2;
	std::valarray<double> f_copy;
	std::valarray<double> tmp;
	std::valarray<double> sol;
	std::map<DomainK, std::valarray<double>> denoms;

	public:
	FftwPatchSolver(DomainCollection &dsc,double lambda=0);
	~FftwPatchSolver();
	void solve(Domain &d, const vector_type &f, vector_type &u, const vector_type &gamma);
	void solve(Domain &d, const vector_type &f, vector_type u, const vector_type &gamma);
	void addDomain(Domain &d);
};
#endif
