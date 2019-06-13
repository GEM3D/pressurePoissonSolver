/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef FFTWPATCHSOLVER_H
#define FFTWPATCHSOLVER_H
#include <Thunderegg/Domain.h>
#include <Thunderegg/PatchSolvers/PatchSolver.h>
#include <Thunderegg/ValVector.h>
#include <bitset>
#include <fftw3.h>
#include <map>
#ifndef DOMAINK
#define DOMAINK
template <size_t D> struct DomainK {
	unsigned long neumann = 0;
	double        h_x     = 0;

	DomainK() {}
	DomainK(const SchurInfo<D> &sinfo)
	{
		this->neumann = sinfo.pinfo->neumann.to_ulong();
		this->h_x     = sinfo.pinfo->spacings[0];
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
	ValVector<D>                                f_copy;
	ValVector<D>                                tmp;
	ValVector<D>                                sol;
	std::array<int, D + 1>                      npow;
	std::map<DomainK<D>, std::valarray<double>> denoms;

	public:
	FftwPatchSolver(Domain<D> &domain, double lambda = 0);
	~FftwPatchSolver();
	void solve(SchurInfo<D> &sinfo, std::shared_ptr<const Vector<D>> f,
	           std::shared_ptr<Vector<D>> u, std::shared_ptr<const Vector<D - 1>> gamma);
	void domainSolve(std::deque<SchurInfo<D>> &domains, std::shared_ptr<const Vector<D>> f,
	                 std::shared_ptr<Vector<D>> u, std::shared_ptr<const Vector<D - 1>> gamma)
	{
		for (SchurInfo<D> &sinfo : domains) {
			solve(sinfo, f, u, gamma);
		}
	}
	void addDomain(SchurInfo<D> &sinfo);
};
template <size_t D> FftwPatchSolver<D>::FftwPatchSolver(Domain<D> &domain, double lambda)
{
	n            = domain.getNs()[0];
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
template <size_t D> void FftwPatchSolver<D>::addDomain(SchurInfo<D> &sinfo)
{
	using namespace std;
	if (!initialized) {
		initialized = true;
		for (size_t i = 0; i <= D; i++) {
			npow[i] = (int) std::pow(n, i);
		}
		std::array<int, D> lengths;
		lengths.fill(n);
		f_copy = ValVector<D>(lengths);
		tmp    = ValVector<D>(lengths);
		sol    = ValVector<D>(lengths);
	}

	int           ns[D];
	fftw_r2r_kind transforms[D];
	fftw_r2r_kind transforms_inv[D];
	if (!plan1.count(sinfo)) {
		for (size_t i = 0; i < D; i++) {
			ns[D - 1 - i] = n;
			// x direction
			if (sinfo.pinfo->isNeumann(2 * i) && sinfo.pinfo->isNeumann(2 * i + 1)) {
				transforms[D - 1 - i]     = FFTW_REDFT10;
				transforms_inv[D - 1 - i] = FFTW_REDFT01;
			} else if (sinfo.pinfo->isNeumann(2 * i)) {
				transforms[D - 1 - i]     = FFTW_REDFT11;
				transforms_inv[D - 1 - i] = FFTW_REDFT11;
			} else if (sinfo.pinfo->isNeumann(2 * i + 1)) {
				transforms[D - 1 - i]     = FFTW_RODFT11;
				transforms_inv[D - 1 - i] = FFTW_RODFT11;
			} else {
				transforms[D - 1 - i]     = FFTW_RODFT10;
				transforms_inv[D - 1 - i] = FFTW_RODFT01;
			}
		}

		plan1[sinfo] = fftw_plan_r2r(D, ns, &f_copy.vec[0], &tmp.vec[0], transforms,
		                             FFTW_MEASURE | FFTW_DESTROY_INPUT);
		plan2[sinfo] = fftw_plan_r2r(D, ns, &tmp.vec[0], &sol.vec[0], transforms_inv,
		                             FFTW_MEASURE | FFTW_DESTROY_INPUT);
	}

	if (!denoms.count(sinfo)) {
		valarray<double> &denom = denoms[sinfo];
		denom.resize(pow(n, D));

		valarray<double> ones(pow(n, D - 1));
		ones = 1;

		for (size_t i = 0; i < D; i++) {
			valarray<size_t> sizes(D - 1);
			sizes = n;
			valarray<size_t> strides(D - 1);
			for (size_t sinfo = 1; sinfo < D; sinfo++) {
				strides[sinfo - 1] = pow(n, (i + sinfo) % D);
			}
			double h = sinfo.pinfo->spacings[i];

			if (sinfo.pinfo->isNeumann(i * 2) && sinfo.pinfo->isNeumann(i * 2 + 1)) {
				for (int xi = 0; xi < n; xi++) {
					denom[gslice(xi * pow(n, i), sizes, strides)]
					-= 4 / (h * h) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
				}
			} else if (sinfo.pinfo->isNeumann(i * 2) || sinfo.pinfo->isNeumann(i * 2 + 1)) {
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
void FftwPatchSolver<D>::solve(SchurInfo<D> &sinfo, std::shared_ptr<const Vector<D>> f,
                               std::shared_ptr<Vector<D>>           u,
                               std::shared_ptr<const Vector<D - 1>> gamma)
{
	using namespace std;
	const LocalData<D> f_view      = f->getLocalData(sinfo.pinfo->local_index);
	LocalData<D>       f_copy_view = f_copy.getLocalData(0);
	LocalData<D>       tmp_view    = tmp.getLocalData(0);

	std::array<int, D> start, end;
	start.fill(0);
	end.fill(n - 1);

	nested_loop<D>(start, end,
	               [&](std::array<int, D> coord) { f_copy_view[coord] = f_view[coord]; });

	for (Side<D> s : Side<D>::getValues()) {
		std::array<int, D - 1> start, end;
		start.fill(0);
		end.fill(n - 1);

		if (sinfo.pinfo->hasNbr(s)) {
			const LocalData<D - 1> gamma_view = gamma->getLocalData(sinfo.getIfaceLocalIndex(s));
			LocalData<D - 1>       slice      = f_copy.getLocalData(0).getSliceOnSide(s);
			double                 h2         = pow(sinfo.pinfo->spacings[s.axis()], 2);
			nested_loop<D - 1>(start, end, [&](std::array<int, D - 1> coord) {
				slice[coord] -= 2.0 / h2 * gamma_view[coord];
			});
		}
	}

	fftw_execute(plan1[sinfo]);

	tmp.vec /= denoms[sinfo];

	if (sinfo.pinfo->neumann.all()) { tmp.vec[0] = 0; }

	fftw_execute(plan2[sinfo]);

	sol.vec /= pow(2.0 * n, D);

	LocalData<D> u_view   = u->getLocalData(sinfo.pinfo->local_index);
	LocalData<D> sol_view = sol.getLocalData(0);
	nested_loop<D>(start, end, [&](std::array<int, D> coord) { u_view[coord] = sol_view[coord]; });
}
extern template class FftwPatchSolver<2>;
extern template class FftwPatchSolver<3>;
#endif
