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

#ifndef DFTPATCHSOLVER_H
#define DFTPATCHSOLVER_H
#include <DomainCollection.h>
#include <PatchSolvers/PatchSolver.h>
#include <Utils.h>
#include <Vector.h>
#include <bitset>
#include <fftw3.h>
#include <map>
#include <valarray>
extern "C" void dgemv_(char &, int &, int &, double &, double *, int &, double *, int &, double &,
                       double *, int &);

enum class DftType { DCT_II, DCT_III, DCT_IV, DST_II, DST_III, DST_IV };
#ifndef DOMAINK
#define DOMAINK
template <size_t D> struct DomainK {
	unsigned long  neumann = 0;
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

template <size_t D> class DftPatchSolver : public PatchSolver<D>
{
	private:
	int                                                                         n;
	bool                                                                        initialized = false;
	static bool                                                                 compareDomains();
	double                                                                      lambda;
	std::map<DomainK<D>, std::array<std::shared_ptr<std::valarray<double>>, D>> plan1;
	std::map<DomainK<D>, std::array<std::shared_ptr<std::valarray<double>>, D>> plan2;
	std::valarray<double>                                                       f_copy;
	std::valarray<double>                                                       tmp;
	std::map<DomainK<D>, std::valarray<double>>                                 eigen_vals;
	std::array<std::shared_ptr<std::valarray<double>>, 6>                       transforms
	= {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}};

	std::array<std::shared_ptr<std::valarray<double>>, D> plan(DftType *types);

	std::shared_ptr<std::valarray<double>> getTransformArray(DftType type);
	void execute_plan(std::array<std::shared_ptr<std::valarray<double>>, D> plan, double *in,
	                  double *out, const bool inverse);

	public:
	DftPatchSolver(DomainCollection<D> &dsc, double lambda = 0);
	void solve(SchurDomain<D> &d, const Vec f, Vec u, const Vec gamma);
	void domainSolve(std::deque<SchurDomain<D>> &domains, const Vec f, Vec u, const Vec gamma)
	{
		for (SchurDomain<D> &d : domains) {
			solve(d, f, u, gamma);
		}
	}
	void addDomain(SchurDomain<D> &d);
};

template <size_t D> inline DftPatchSolver<D>::DftPatchSolver(DomainCollection<D> &dc, double lambda)
{
	n            = dc.getN();
	this->lambda = lambda;
}
template <size_t D> inline void DftPatchSolver<D>::addDomain(SchurDomain<D> &d)
{
	using namespace std;
	if (!initialized) {
		initialized = true;
		f_copy.resize(pow(n, D));
		tmp.resize(pow(n, D));
	}

	DftType transforms[D];
	DftType transforms_inv[D];
	if (!plan1.count(d)) {
		for (size_t i = 0; i < D; i++) {
			// x direction
			if (d.isNeumann(2 * i) && d.isNeumann(2 * i + 1)) {
				transforms[i]     = DftType::DCT_II;
				transforms_inv[i] = DftType::DCT_III;
			} else if (d.isNeumann(2 * i)) {
				transforms[i]     = DftType::DCT_IV;
				transforms_inv[i] = DftType::DCT_IV;
			} else if (d.isNeumann(2 * i + 1)) {
				transforms[i]     = DftType::DST_IV;
				transforms_inv[i] = DftType::DST_IV;
			} else {
				transforms[i]     = DftType::DST_II;
				transforms_inv[i] = DftType::DST_III;
			}
		}

		plan1[d] = plan(transforms);
		plan2[d] = plan(transforms_inv);
	}

	if (!eigen_vals.count(d)) {
		valarray<double> &denom = eigen_vals[d];
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
inline void DftPatchSolver<D>::solve(SchurDomain<D> &d, const Vec f, Vec u, const Vec gamma)
{
	using namespace std;
	using namespace Utils;
	const double *f_view, *gamma_view;
	VecGetArrayRead(f, &f_view);
	VecGetArrayRead(gamma, &gamma_view);

	int start = d.local_index * pow(n, D);
	for (int i = 0; i < (const int) pow(n, D); i++) {
		f_copy[i] = f_view[start + i];
	}

	for (Side<D> s : Side<D>::getValues()) {
		if (d.hasNbr(s)) {
			int          idx = pow(n, D - 1) * d.getIfaceLocalIndex(s);
			Slice<D - 1> sl  = getSlice<D - 1>(&f_copy[0], n, s);
			double       h2  = pow(d.domain.lengths[s.toInt() / 2] / n, 2);
			int          strides[D - 1];
			for (size_t i = 0; i < D - 1; i++) {
				strides[i] = pow(n, i);
			}
			for (int i = 0; i < (const int) pow(n, D - 1); i++) {
				std::array<int, D - 1> coord;
				for (size_t x = 0; x < D - 1; x++) {
					coord[x] = (i / strides[x]) % n;
				}
				sl(coord) -= 2.0 / h2 * gamma_view[idx + i];
			}
		}
	}

	execute_plan(plan1[d], &f_copy[0], &tmp[0], false);

	tmp /= eigen_vals[d];

	if (d.neumann.all()) { tmp[0] = 0; }

	double *u_view;
	VecGetArray(u, &u_view);
	double *u_local_view = u_view + start;
	execute_plan(plan2[d], &tmp[0], u_local_view, true);

	double scale = pow(2.0 / n, D);
	for (int i = 0; i < (const int) pow(n, D); i++) {
		u_local_view[i] *= scale;
	}
	VecRestoreArray(u, &u_view);
	VecRestoreArrayRead(f, &f_view);
	VecRestoreArrayRead(gamma, &gamma_view);
}
template <size_t D>
inline std::array<std::shared_ptr<std::valarray<double>>, D> DftPatchSolver<D>::plan(DftType *types)
{
	std::array<std::shared_ptr<std::valarray<double>>, D> retval;
	for (size_t i = 0; i < D; i++) {
		retval[i] = getTransformArray(types[i]);
	}
	return retval;
}
template <size_t D>
inline std::shared_ptr<std::valarray<double>> DftPatchSolver<D>::getTransformArray(DftType type)
{
	using namespace std;
	int idx = static_cast<int>(type);

	shared_ptr<valarray<double>> matrix_ptr;

	if (transforms[idx] == nullptr) {
		matrix_ptr.reset(new valarray<double>(n * n));
		valarray<double> &matrix = *matrix_ptr;
		switch (type) {
			case DftType::DCT_II:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * (i * (j + 0.5)));
					}
				}
				break;
			case DftType::DCT_III:
				for (int i = 0; i < n; i++) {
					matrix[i * n] = 0.5;
				}
				for (int j = 1; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * ((i + 0.5) * j));
					}
				}
				break;
			case DftType::DCT_IV:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * ((i + 0.5) * (j + 0.5)));
					}
				}
				break;
			case DftType::DST_II:
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 1) * (j + 0.5)));
					}
				}
				break;
			case DftType::DST_III:
				for (int i = 0; i < n; i += 2) {
					matrix[i * n + n - 1] = 0.5;
				}
				for (int i = 1; i < n; i += 2) {
					matrix[i * n + n - 1] = -0.5;
				}
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n - 1; j++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 0.5) * (j + 1)));
					}
				}
				break;
			case DftType::DST_IV:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 0.5) * (j + 0.5)));
					}
				}
				break;
		}
	} else {
		matrix_ptr = transforms[idx];
	}
	return matrix_ptr;
}
template <size_t D>
inline void
DftPatchSolver<D>::execute_plan(std::array<std::shared_ptr<std::valarray<double>>, D> plan,
                                double *in, double *out, const bool inverse)
{
	double *prev_result = in;
	int     strides[D];
	for (size_t i = 0; i < D; i++) {
		strides[i] = pow(n, i);
	}
	for (size_t dim = 0; dim < D; dim++) {
		int other_strides[D - 1];
		for (size_t i = 0; i < dim; i++) {
			other_strides[i] = strides[i];
		}
		int dft_stride = strides[dim];
		for (size_t i = dim + 1; i < D; i++) {
			other_strides[i - 1] = strides[i];
		}

		std::valarray<double> &matrix = *plan[dim];
		double *               new_result;
		if (dim != D-1) {
			new_result = new double[(int)pow(n,D)];
		} else {
			new_result = out;
		}

		for (int i = 0; i < (const int) pow(n, D - 1); i++) {
			std::array<int, D - 1> coord;
			for (size_t x = 0; x < D - 1; x++) {
				coord[x] = (i / strides[x]) % n;
			}
			int    idx = std::inner_product(other_strides, other_strides + D - 1, coord.begin(), 0);
			char   T   = 'T';
			double one = 1;
			double zero = 0;
			dgemv_(T, n, n, one, &matrix[0], n, &prev_result[idx], dft_stride, zero,
			       &new_result[idx], dft_stride);
		}
		if (dim != 0) { delete[] prev_result; }
		prev_result = new_result;
	}
}
#endif
