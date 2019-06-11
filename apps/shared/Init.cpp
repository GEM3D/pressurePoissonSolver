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

#include "Init.h"
#include <algorithm>
using namespace std;
void getXYZ(const PatchInfo<3> &d, const int &xi, const int &yi, const int &zi, double &x,
            double &y, double &z)
{
	double h_x = d.spacings[0];
	double h_y = d.spacings[1];
	double h_z = d.spacings[2];
	if (xi == -1) {
		x = d.starts[0];
	} else if (xi == d.ns[0]) {
		x = d.starts[0] + d.spacings[0] * d.ns[0];
	} else {
		x = d.starts[0] + h_x / 2.0 + d.spacings[0] * xi;
	}
	if (yi == -1) {
		y = d.starts[1];
	} else if (yi == d.ns[1]) {
		y = d.starts[1] + d.spacings[1] * d.ns[1];
	} else {
		y = d.starts[1] + h_y / 2.0 + d.spacings[1] * yi;
	}
	if (zi == -1) {
		z = d.starts[2];
	} else if (zi == d.ns[2]) {
		z = d.starts[2] + d.spacings[2] * d.ns[2];
	} else {
		z = d.starts[2] + h_z / 2.0 + d.spacings[2] * zi;
	}
}
inline int index(PatchInfo<3> &d, const int &xi, const int &yi, const int &zi)
{
	return xi + yi * d.ns[0] + zi * d.ns[0] * d.ns[1];
}
void Init::initNeumann(DomainCollection<3> &dc, Vec f, Vec exact,
                       function<double(double, double, double)> ffun,
                       function<double(double, double, double)> efun,
                       function<double(double, double, double)> nfunx,
                       function<double(double, double, double)> nfuny,
                       function<double(double, double, double)> nfunz)
{
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &p : dc.domains) {
		PatchInfo<3> &d = *p.second;

		double *f_vals     = f_ptr + d.local_index * d.ns[0] * d.ns[1] * d.ns[2];
		double *exact_vals = exact_ptr + d.local_index * d.ns[0] * d.ns[1] * d.ns[2];

		// Generate RHS vector
		for (int zi = 0; zi < d.ns[2]; zi++) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, zi, x, y, z);
					f_vals[index(d, xi, yi, zi)]     = ffun(x, y, z);
					exact_vals[index(d, xi, yi, zi)] = efun(x, y, z);
				}
			}
		}
		// apply boundaries
		double h_x = d.spacings[0];
		double h_y = d.spacings[1];
		double h_z = d.spacings[2];
		// west
		if (!d.hasNbr(Side<3>::west)) {
			for (int zi = 0; zi < d.ns[2]; zi++) {
				for (int yi = 0; yi < d.ns[1]; yi++) {
					double x, y, z;
					getXYZ(d, -1, yi, zi, x, y, z);
					f_vals[index(d, 0, yi, zi)] += nfunx(x, y, z) / h_x;
				}
			}
		}
		// east
		if (!d.hasNbr(Side<3>::east)) {
			for (int zi = 0; zi < d.ns[2]; zi++) {
				for (int yi = 0; yi < d.ns[1]; yi++) {
					double x, y, z;
					getXYZ(d, d.ns[0], yi, zi, x, y, z);
					f_vals[index(d, d.ns[0] - 1, yi, zi)] -= nfunx(x, y, z) / h_x;
				}
			}
		}
		// south
		if (!d.hasNbr(Side<3>::south)) {
			for (int zi = 0; zi < d.ns[2]; zi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, -1, zi, x, y, z);
					f_vals[index(d, xi, 0, zi)] += nfuny(x, y, z) / h_y;
				}
			}
		}
		// north
		if (!d.hasNbr(Side<3>::north)) {
			for (int zi = 0; zi < d.ns[2]; zi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, d.ns[1], zi, x, y, z);
					f_vals[index(d, xi, d.ns[1] - 1, zi)] -= nfuny(x, y, z) / h_y;
				}
			}
		}
		// bottom
		if (!d.hasNbr(Side<3>::bottom)) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, -1, x, y, z);
					f_vals[index(d, xi, yi, 0)] += nfunz(x, y, z) / h_z;
				}
			}
		}
		// top
		if (!d.hasNbr(Side<3>::top)) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, d.ns[2], x, y, z);
					f_vals[index(d, xi, yi, d.ns[2] - 1)] -= nfunz(x, y, z) / h_z;
				}
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::initDirichlet(DomainCollection<3> &dc, Vec f, Vec exact,
                         function<double(double, double, double)> ffun,
                         function<double(double, double, double)> efun)
{
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &p : dc.domains) {
		PatchInfo<3> &d = *p.second;

		double *f_vals     = f_ptr + d.local_index * d.ns[0] * d.ns[1] * d.ns[2];
		double *exact_vals = exact_ptr + d.local_index * d.ns[0] * d.ns[1] * d.ns[2];
		// Generate RHS vector
		double h_x = d.spacings[0];
		double h_y = d.spacings[1];
		double h_z = d.spacings[2];
		h_x *= h_x;
		h_y *= h_y;
		h_z *= h_z;
		for (int zi = 0; zi < d.ns[2]; zi++) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, zi, x, y, z);
					f_vals[index(d, xi, yi, zi)]     = ffun(x, y, z);
					exact_vals[index(d, xi, yi, zi)] = efun(x, y, z);
				}
			}
		}
		// apply boundaries
		// west
		if (!d.hasNbr(Side<3>::west)) {
			for (int zi = 0; zi < d.ns[2]; zi++) {
				for (int yi = 0; yi < d.ns[1]; yi++) {
					double x, y, z;
					getXYZ(d, -1, yi, zi, x, y, z);
					f_vals[index(d, 0, yi, zi)] -= 2 * efun(x, y, z) / h_x;
				}
			}
		}
		// east
		if (!d.hasNbr(Side<3>::east)) {
			for (int zi = 0; zi < d.ns[2]; zi++) {
				for (int yi = 0; yi < d.ns[1]; yi++) {
					double x, y, z;
					getXYZ(d, d.ns[0], yi, zi, x, y, z);
					f_vals[index(d, d.ns[0] - 1, yi, zi)] -= 2 * efun(x, y, z) / h_x;
				}
			}
		}
		// south
		if (!d.hasNbr(Side<3>::south)) {
			for (int zi = 0; zi < d.ns[2]; zi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, -1, zi, x, y, z);
					f_vals[index(d, xi, 0, zi)] -= 2 * efun(x, y, z) / h_y;
				}
			}
		}
		// north
		if (!d.hasNbr(Side<3>::north)) {
			for (int zi = 0; zi < d.ns[2]; zi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, d.ns[1], zi, x, y, z);
					f_vals[index(d, xi, d.ns[1] - 1, zi)] -= 2 * efun(x, y, z) / h_y;
				}
			}
		}
		// bottom
		if (!d.hasNbr(Side<3>::bottom)) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, -1, x, y, z);
					f_vals[index(d, xi, yi, 0)] -= 2 * efun(x, y, z) / h_z;
				}
			}
		}
		// top
		if (!d.hasNbr(Side<3>::top)) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				for (int xi = 0; xi < d.ns[0]; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, d.ns[2], x, y, z);
					f_vals[index(d, xi, yi, d.ns[2] - 1)] -= 2 * efun(x, y, z) / h_z;
				}
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::initNeumann2d(DomainCollection<2> &dc, Vec f, Vec exact,
                         function<double(double, double)> ffun,
                         function<double(double, double)> efun,
                         function<double(double, double)> nfunx,
                         function<double(double, double)> nfuny)
{
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &p : dc.domains) {
		PatchInfo<2> &d = *p.second;

		double *f_vals     = f_ptr + d.local_index * d.ns[0] * d.ns[1];
		double *exact_vals = exact_ptr + d.local_index * d.ns[0] * d.ns[1];

		// Generate RHS vector
		double h_x = d.spacings[0];
		double h_y = d.spacings[1];
		for (int yi = 0; yi < d.ns[1]; yi++) {
			for (int xi = 0; xi < d.ns[0]; xi++) {
				double x                      = d.starts[0] + h_x / 2.0 + h_x * xi;
				double y                      = d.starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * d.ns[0] + xi]     = ffun(x, y);
				exact_vals[yi * d.ns[0] + xi] = efun(x, y);
			}
		}
		// apply boundaries
		// west
		if (!d.hasNbr(Side<2>::west)) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				double y = d.starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * d.ns[0]] += nfunx(d.starts[0], y) / h_x;
			}
		}
		// east
		if (!d.hasNbr(Side<2>::east)) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				double y = d.starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * d.ns[0] + d.ns[0] - 1] -= nfunx(d.starts[0] + h_x * d.ns[0], y) / h_x;
			}
		}
		// south
		if (!d.hasNbr(Side<2>::south)) {
			for (int xi = 0; xi < d.ns[0]; xi++) {
				double x = d.starts[0] + h_x / 2.0 + h_x * xi;
				f_vals[xi] += nfuny(x, d.starts[1]) / h_y;
			}
		}
		// north
		if (!d.hasNbr(Side<2>::north)) {
			for (int xi = 0; xi < d.ns[0]; xi++) {
				double x = d.starts[0] + h_x / 2.0 + h_x * xi;
				f_vals[d.ns[0] * (d.ns[1] - 1) + xi] -= nfuny(x, d.starts[1] + h_y * d.ns[1]) / h_y;
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::initDirichlet2d(DomainCollection<2> &dc, Vec f, Vec exact,
                           function<double(double, double)> ffun,
                           function<double(double, double)> efun)
{
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &p : dc.domains) {
		PatchInfo<2> &d = *p.second;

		double *f_vals     = f_ptr + d.local_index * d.ns[0] * d.ns[1];
		double *exact_vals = exact_ptr + d.local_index * d.ns[0] * d.ns[1];
		// Generate RHS vector
		double h_x = d.spacings[0];
		double h_y = d.spacings[1];
		for (int yi = 0; yi < d.ns[1]; yi++) {
			for (int xi = 0; xi < d.ns[0]; xi++) {
				double x                      = d.starts[0] + h_x / 2.0 + h_x * xi;
				double y                      = d.starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * d.ns[0] + xi]     = ffun(x, y);
				exact_vals[yi * d.ns[0] + xi] = efun(x, y);
			}
		}
		// apply boundaries
		// west
		if (!d.hasNbr(Side<2>::west)) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				double y = d.starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * d.ns[0]] -= efun(d.starts[0], y) * 2 / (h_x * h_x);
			}
		}
		// east
		if (!d.hasNbr(Side<2>::east)) {
			for (int yi = 0; yi < d.ns[1]; yi++) {
				double y = d.starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * d.ns[0] + d.ns[0] - 1]
				-= efun(d.starts[0] + h_x * d.ns[0], y) * 2 / (h_x * h_x);
			}
		}
		// south
		if (!d.hasNbr(Side<2>::south)) {
			for (int xi = 0; xi < d.ns[0]; xi++) {
				double x = d.starts[0] + h_x / 2.0 + h_x * xi;
				f_vals[xi] -= efun(x, d.starts[1]) * 2 / (h_y * h_y);
			}
		}
		// north
		if (!d.hasNbr(Side<2>::north)) {
			for (int xi = 0; xi < d.ns[0]; xi++) {
				double x = d.starts[0] + h_x / 2.0 + h_x * xi;
				f_vals[d.ns[0] * (d.ns[1] - 1) + xi]
				-= efun(x, d.starts[1] + h_y * d.ns[1]) * 2 / (h_y * h_y);
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::fillSolution2d(DomainCollection<2> &dc, Vec u,
                          function<double(double, double, double)> fun, double time)
{
	double *vec;
	VecGetArray(u, &vec);
	for (auto &p : dc.domains) {
		PatchInfo<2> &d = *p.second;

		double *f = vec + d.local_index * d.ns[0] * d.ns[1];

		double h_x = d.spacings[0];
		double h_y = d.spacings[1];
		for (int yi = 0; yi < d.ns[1]; yi++) {
			for (int xi = 0; xi < d.ns[0]; xi++) {
				double x             = d.starts[0] + h_x / 2.0 + h_x * xi;
				double y             = d.starts[1] + h_y / 2.0 + h_y * yi;
				f[yi * d.ns[0] + xi] = fun(x, y, time);
			}
		}
	}
	VecRestoreArray(u, &vec);
}
