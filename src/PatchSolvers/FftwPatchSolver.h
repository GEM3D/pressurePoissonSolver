#ifndef FFTWPATCHSOLVER_H
#define FFTWPATCHSOLVER_H
#include "DomainCollection.h"
#include "PatchSolvers/PatchSolver.h"
#include "Utils.h"
#include <bitset>
#include <fftw3.h>
#include <map>
#include <valarray>
#ifndef DOMAINK
#define DOMAINK
template <size_t D> struct DomainK {
	ulong  neumann = 0;
	double h_x     = 0;

	DomainK() {}
	DomainK(const SchurDomain<D> &d)
	{
		this->neumann = d.neumann.to_ulong();
		this->h_x     = d.domain.lengths[0];
	}
	friend bool operator<(const DomainK &l, const DomainK &r)
	{
		return std::tie(l.neumann, l.h_x) < std::tie(r.neumann, r.h_x);
	}
};
#endif

template <size_t D> class FftwPatchSolver : public PatchSolver<D>
{
	private:
	int                                         n;
	bool                                        initialized = false;
	static bool                                 compareDomains();
	double                                      lambda;
	std::map<DomainK<D>, fftw_plan>             plan1;
	std::map<DomainK<D>, fftw_plan>             plan2;
	std::valarray<double>                       f_copy;
	std::valarray<double>                       tmp;
	std::valarray<double>                       sol;
	std::map<DomainK<D>, std::valarray<double>> denoms;

	public:
	FftwPatchSolver(DomainCollection<D> &dsc, double lambda = 0);
	~FftwPatchSolver();
	void solve(SchurDomain<D> &d, const Vec f, Vec u, const Vec gamma);
	void domainSolve(std::deque<SchurDomain<D>> &domains, const Vec f, Vec u, const Vec gamma)
	{
		for (SchurDomain<D> &d : domains) {
			solve(d, f, u, gamma);
		}
	}
	void addDomain(SchurDomain<D> &d);
};
template <size_t D> FftwPatchSolver<D>::FftwPatchSolver(DomainCollection<D> &dc, double lambda)
{
	n            = dc.getN();
	this->lambda = lambda;
}
template <size_t D> FftwPatchSolver<D>::~FftwPatchSolver()
{
	for (auto p : plan1) {
		fftw_destroy_plan(p.second);
	}
	for (auto p : plan2) {
		fftw_destroy_plan(p.second);
	}
}
template <size_t D> void FftwPatchSolver<D>::addDomain(SchurDomain<D> &d)
{
	using namespace std;
	if (!initialized) {
		initialized = true;
		f_copy.resize(pow(n, D));
		tmp.resize(pow(n, D));
		sol.resize(pow(n, D));
	}

	int           ns[D];
	fftw_r2r_kind transforms[D];
	fftw_r2r_kind transforms_inv[D];
	if (!plan1.count(d)) {
		for (size_t i = 0; i < D; i++) {
			ns[D-1-i] = n;
			// x direction
			if (d.isNeumann(2 * i) && d.isNeumann(2 * i + 1)) {
				transforms[D-1-i]     = FFTW_REDFT10;
				transforms_inv[D-1-i] = FFTW_REDFT01;
			} else if (d.isNeumann(2 * i)) {
				transforms[D-1-i]     = FFTW_REDFT11;
				transforms_inv[D-1-i] = FFTW_REDFT11;
			} else if (d.isNeumann(2 * i + 1)) {
				transforms[D-1-i]     = FFTW_RODFT11;
				transforms_inv[D-1-i] = FFTW_RODFT11;
			} else {
				transforms[D-1-i]     = FFTW_RODFT10;
				transforms_inv[D-1-i] = FFTW_RODFT01;
			}
		}

		plan1[d]
		= fftw_plan_r2r(D, ns, &f_copy[0], &tmp[0], transforms, FFTW_MEASURE | FFTW_DESTROY_INPUT);
		plan2[d]
		= fftw_plan_r2r(D, ns, &tmp[0], &sol[0], transforms_inv, FFTW_MEASURE | FFTW_DESTROY_INPUT);
	}

	if (!denoms.count(d)) {
		valarray<double> &denom = denoms[d];
		denom.resize(pow(n, D));

		valarray<double> ones(pow(n, D - 1));
		ones = 1;

		for (size_t i = 0; i < D; i++) {
			valarray<size_t> sizes(D - 1);
			sizes = n;
			valarray<size_t> strides(D - 1);
			for (size_t d = 1; d < D; d++) {
				strides[d - 1] = pow(n, (i + d) % D);
			}
			double h = d.domain.lengths[i] / n;

			if (d.isNeumann(i * 2) && d.isNeumann(i * 2 + 1)) {
				for (int xi = 0; xi < n; xi++) {
					denom[gslice(xi * pow(n, i), sizes, strides)]
					-= 4 / (h * h) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
				}
			} else if (d.isNeumann(i * 2) || d.isNeumann(i * 2 + 1)) {
				for (int xi = 0; xi < n; xi++) {
					denom[gslice(xi * pow(n, i), sizes, strides)]
					-= 4 / (h * h) * pow(sin((xi + 0.5) * M_PI / (2 * n)), 2) * ones;
				}
			} else {
				for (int xi = 0; xi < n; xi++) {
					denom[gslice(xi * pow(n, i), sizes, strides)]
					-= 4 / (h * h) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
				}
			}
		}

		denom += lambda;
	}
}
template <size_t D>
void FftwPatchSolver<D>::solve(SchurDomain<D> &d, const Vec f, Vec u, const Vec gamma)
{
	using namespace std;
	using namespace Utils;
	const double *f_view, *gamma_view;
	VecGetArrayRead(f, &f_view);
	VecGetArrayRead(gamma, &gamma_view);

	int start = d.local_index * pow(n, D);
	for (int i = 0; i < pow(n, D); i++) {
		f_copy[i] = f_view[start + i];
	}

	for (Side<D> s : Side<D>::getValues()) {
		if (d.hasNbr(s)) {
			int          idx = pow(n, D - 1) * d.getIfaceLocalIndex(s);
			Slice<D - 1> sl  = getSlice<D - 1>(&f_copy[0], n, s);
			double       h2  = pow(d.domain.lengths[s.toInt() / 2] / n, 2);
			for (int i = 0; i < (int) pow(n, D - 1); i++) {
				std::array<int, D - 1> coord;
				for (size_t x = 0; x < D - 1; x++) {
					coord[x] = (i / (int) pow(n, x)) % n;
				}
				sl(coord) -= 2.0 / h2 * gamma_view[idx + i];
			}
		}
	}

	fftw_execute(plan1[d]);

	tmp /= denoms[d];

	if (d.neumann.all()) { tmp[0] = 0; }

	fftw_execute(plan2[d]);

	sol /= pow(2.0 * n, D);

	double *u_view;
	VecGetArray(u, &u_view);
	for (int i = 0; i < pow(n, D); i++) {
		u_view[start + i] = sol[i];
	}
	VecRestoreArray(u, &u_view);
	VecRestoreArrayRead(f, &f_view);
	VecRestoreArrayRead(gamma, &gamma_view);
}
#endif
