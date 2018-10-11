#ifndef DFTPATCHSOLVER_H
#define DFTPATCHSOLVER_H
#include "DomainCollection.h"
#include "PatchSolvers/PatchSolver.h"
#include <bitset>
#include <fftw3.h>
#include <map>
#include <valarray>
enum class DftType { DCT_II, DCT_III, DCT_IV, DST_II, DST_III, DST_IV };
#ifndef DOMAINK
#define DOMAINK
struct DomainK {
	ulong  neumann = 0;
	double h_x     = 0;
	double h_y     = 0;

	DomainK() {}
	DomainK(const SchurDomain &d)
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
#endif

class DftPatchSolver : public PatchSolver
{
	private:
	int                                                                      n;
	bool                                                                     initialized = false;
	static bool                                                              compareDomains();
	double                                                                   lambda;
	std::map<DomainK, std::array<std::shared_ptr<std::valarray<double>>, 3>> plan1;
	std::map<DomainK, std::array<std::shared_ptr<std::valarray<double>>, 3>> plan2;
	std::valarray<double>                                                    f_copy;
	std::valarray<double>                                                    tmp;
	std::map<DomainK, std::valarray<double>>                                 eigen_vals;
	std::array<std::shared_ptr<std::valarray<double>>, 6>                    transforms
	= {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}};

	std::array<std::shared_ptr<std::valarray<double>>, 3> plan(DftType x_type, DftType y_type,
	                                                           DftType z_type);

	std::shared_ptr<std::valarray<double>> getTransformArray(DftType type);
	void execute_plan(std::array<std::shared_ptr<std::valarray<double>>, 3> plan, double *in,
	                  double *out, const bool inverse);

	public:
	DftPatchSolver(DomainCollection &dsc, double lambda = 0);
	void solve(SchurDomain &d, const Vec f, Vec u, const Vec gamma);
	void domainSolve(std::deque<SchurDomain> &domains, const Vec f, Vec u, const Vec gamma){
        for(SchurDomain &d:domains){
            solve(d,f,u,gamma);
        }
    }
	void addDomain(SchurDomain &d);
};
#endif