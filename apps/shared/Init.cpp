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
static void getXYZ(shared_ptr<PatchInfo<3>> pinfo, int xi, int yi, int zi, double &x, double &y,
                   double &z)
{
	double h_x = pinfo->spacings[0];
	double h_y = pinfo->spacings[1];
	double h_z = pinfo->spacings[2];
	if (xi == -1) {
		x = pinfo->starts[0];
	} else if (xi == pinfo->ns[0]) {
		x = pinfo->starts[0] + pinfo->spacings[0] * pinfo->ns[0];
	} else {
		x = pinfo->starts[0] + h_x / 2.0 + pinfo->spacings[0] * xi;
	}
	if (yi == -1) {
		y = pinfo->starts[1];
	} else if (yi == pinfo->ns[1]) {
		y = pinfo->starts[1] + pinfo->spacings[1] * pinfo->ns[1];
	} else {
		y = pinfo->starts[1] + h_y / 2.0 + pinfo->spacings[1] * yi;
	}
	if (zi == -1) {
		z = pinfo->starts[2];
	} else if (zi == pinfo->ns[2]) {
		z = pinfo->starts[2] + pinfo->spacings[2] * pinfo->ns[2];
	} else {
		z = pinfo->starts[2] + h_z / 2.0 + pinfo->spacings[2] * zi;
	}
}
inline int index(shared_ptr<PatchInfo<3>> pinfo, int xi, int yi, int zi)
{
	return xi + yi * pinfo->ns[0] + zi * pinfo->ns[0] * pinfo->ns[1];
}
void Init::initNeumann(Domain<3> &domain, Vec f, Vec exact,
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
	for (auto &pinfo : domain.getPatchInfoVector()) {
		double *f_vals = f_ptr + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1] * pinfo->ns[2];
		double *exact_vals
		= exact_ptr + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1] * pinfo->ns[2];

		// Generate RHS vector
		for (int zi = 0; zi < pinfo->ns[2]; zi++) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, yi, zi, x, y, z);
					f_vals[index(pinfo, xi, yi, zi)]     = ffun(x, y, z);
					exact_vals[index(pinfo, xi, yi, zi)] = efun(x, y, z);
				}
			}
		}
		// apply boundaries
		double h_x = pinfo->spacings[0];
		double h_y = pinfo->spacings[1];
		double h_z = pinfo->spacings[2];
		// west
		if (!pinfo->hasNbr(Side<3>::west)) {
			for (int zi = 0; zi < pinfo->ns[2]; zi++) {
				for (int yi = 0; yi < pinfo->ns[1]; yi++) {
					double x, y, z;
					getXYZ(pinfo, -1, yi, zi, x, y, z);
					f_vals[index(pinfo, 0, yi, zi)] += nfunx(x, y, z) / h_x;
				}
			}
		}
		// east
		if (!pinfo->hasNbr(Side<3>::east)) {
			for (int zi = 0; zi < pinfo->ns[2]; zi++) {
				for (int yi = 0; yi < pinfo->ns[1]; yi++) {
					double x, y, z;
					getXYZ(pinfo, pinfo->ns[0], yi, zi, x, y, z);
					f_vals[index(pinfo, pinfo->ns[0] - 1, yi, zi)] -= nfunx(x, y, z) / h_x;
				}
			}
		}
		// south
		if (!pinfo->hasNbr(Side<3>::south)) {
			for (int zi = 0; zi < pinfo->ns[2]; zi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, -1, zi, x, y, z);
					f_vals[index(pinfo, xi, 0, zi)] += nfuny(x, y, z) / h_y;
				}
			}
		}
		// north
		if (!pinfo->hasNbr(Side<3>::north)) {
			for (int zi = 0; zi < pinfo->ns[2]; zi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, pinfo->ns[1], zi, x, y, z);
					f_vals[index(pinfo, xi, pinfo->ns[1] - 1, zi)] -= nfuny(x, y, z) / h_y;
				}
			}
		}
		// bottom
		if (!pinfo->hasNbr(Side<3>::bottom)) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, yi, -1, x, y, z);
					f_vals[index(pinfo, xi, yi, 0)] += nfunz(x, y, z) / h_z;
				}
			}
		}
		// top
		if (!pinfo->hasNbr(Side<3>::top)) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, yi, pinfo->ns[2], x, y, z);
					f_vals[index(pinfo, xi, yi, pinfo->ns[2] - 1)] -= nfunz(x, y, z) / h_z;
				}
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::initDirichlet(Domain<3> &domain, Vec f, Vec exact,
                         function<double(double, double, double)> ffun,
                         function<double(double, double, double)> efun)
{
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &pinfo : domain.getPatchInfoVector()) {
		double *f_vals = f_ptr + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1] * pinfo->ns[2];
		double *exact_vals
		= exact_ptr + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1] * pinfo->ns[2];
		// Generate RHS vector
		double h_x = pinfo->spacings[0];
		double h_y = pinfo->spacings[1];
		double h_z = pinfo->spacings[2];
		h_x *= h_x;
		h_y *= h_y;
		h_z *= h_z;
		for (int zi = 0; zi < pinfo->ns[2]; zi++) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, yi, zi, x, y, z);
					f_vals[index(pinfo, xi, yi, zi)]     = ffun(x, y, z);
					exact_vals[index(pinfo, xi, yi, zi)] = efun(x, y, z);
				}
			}
		}
		// apply boundaries
		// west
		if (!pinfo->hasNbr(Side<3>::west)) {
			for (int zi = 0; zi < pinfo->ns[2]; zi++) {
				for (int yi = 0; yi < pinfo->ns[1]; yi++) {
					double x, y, z;
					getXYZ(pinfo, -1, yi, zi, x, y, z);
					f_vals[index(pinfo, 0, yi, zi)] -= 2 * efun(x, y, z) / h_x;
				}
			}
		}
		// east
		if (!pinfo->hasNbr(Side<3>::east)) {
			for (int zi = 0; zi < pinfo->ns[2]; zi++) {
				for (int yi = 0; yi < pinfo->ns[1]; yi++) {
					double x, y, z;
					getXYZ(pinfo, pinfo->ns[0], yi, zi, x, y, z);
					f_vals[index(pinfo, pinfo->ns[0] - 1, yi, zi)] -= 2 * efun(x, y, z) / h_x;
				}
			}
		}
		// south
		if (!pinfo->hasNbr(Side<3>::south)) {
			for (int zi = 0; zi < pinfo->ns[2]; zi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, -1, zi, x, y, z);
					f_vals[index(pinfo, xi, 0, zi)] -= 2 * efun(x, y, z) / h_y;
				}
			}
		}
		// north
		if (!pinfo->hasNbr(Side<3>::north)) {
			for (int zi = 0; zi < pinfo->ns[2]; zi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, pinfo->ns[1], zi, x, y, z);
					f_vals[index(pinfo, xi, pinfo->ns[1] - 1, zi)] -= 2 * efun(x, y, z) / h_y;
				}
			}
		}
		// bottom
		if (!pinfo->hasNbr(Side<3>::bottom)) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, yi, -1, x, y, z);
					f_vals[index(pinfo, xi, yi, 0)] -= 2 * efun(x, y, z) / h_z;
				}
			}
		}
		// top
		if (!pinfo->hasNbr(Side<3>::top)) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				for (int xi = 0; xi < pinfo->ns[0]; xi++) {
					double x, y, z;
					getXYZ(pinfo, xi, yi, pinfo->ns[2], x, y, z);
					f_vals[index(pinfo, xi, yi, pinfo->ns[2] - 1)] -= 2 * efun(x, y, z) / h_z;
				}
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::initNeumann2d(Domain<2> &domain, Vec f, Vec exact, function<double(double, double)> ffun,
                         function<double(double, double)> efun,
                         function<double(double, double)> nfunx,
                         function<double(double, double)> nfuny)
{
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &pinfo : domain.getPatchInfoVector()) {
		double *f_vals     = f_ptr + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1];
		double *exact_vals = exact_ptr + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1];

		// Generate RHS vector
		double h_x = pinfo->spacings[0];
		double h_y = pinfo->spacings[1];
		for (int yi = 0; yi < pinfo->ns[1]; yi++) {
			for (int xi = 0; xi < pinfo->ns[0]; xi++) {
				double x                           = pinfo->starts[0] + h_x / 2.0 + h_x * xi;
				double y                           = pinfo->starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * pinfo->ns[0] + xi]     = ffun(x, y);
				exact_vals[yi * pinfo->ns[0] + xi] = efun(x, y);
			}
		}
		// apply boundaries
		// west
		if (!pinfo->hasNbr(Side<2>::west)) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				double y = pinfo->starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * pinfo->ns[0]] += nfunx(pinfo->starts[0], y) / h_x;
			}
		}
		// east
		if (!pinfo->hasNbr(Side<2>::east)) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				double y = pinfo->starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * pinfo->ns[0] + pinfo->ns[0] - 1]
				-= nfunx(pinfo->starts[0] + h_x * pinfo->ns[0], y) / h_x;
			}
		}
		// south
		if (!pinfo->hasNbr(Side<2>::south)) {
			for (int xi = 0; xi < pinfo->ns[0]; xi++) {
				double x = pinfo->starts[0] + h_x / 2.0 + h_x * xi;
				f_vals[xi] += nfuny(x, pinfo->starts[1]) / h_y;
			}
		}
		// north
		if (!pinfo->hasNbr(Side<2>::north)) {
			for (int xi = 0; xi < pinfo->ns[0]; xi++) {
				double x = pinfo->starts[0] + h_x / 2.0 + h_x * xi;
				f_vals[pinfo->ns[0] * (pinfo->ns[1] - 1) + xi]
				-= nfuny(x, pinfo->starts[1] + h_y * pinfo->ns[1]) / h_y;
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::initDirichlet2d(Domain<2> &domain, Vec f, Vec exact,
                           function<double(double, double)> ffun,
                           function<double(double, double)> efun)
{
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &pinfo : domain.getPatchInfoVector()) {
		double *f_vals     = f_ptr + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1];
		double *exact_vals = exact_ptr + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1];
		// Generate RHS vector
		double h_x = pinfo->spacings[0];
		double h_y = pinfo->spacings[1];
		for (int yi = 0; yi < pinfo->ns[1]; yi++) {
			for (int xi = 0; xi < pinfo->ns[0]; xi++) {
				double x                           = pinfo->starts[0] + h_x / 2.0 + h_x * xi;
				double y                           = pinfo->starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * pinfo->ns[0] + xi]     = ffun(x, y);
				exact_vals[yi * pinfo->ns[0] + xi] = efun(x, y);
			}
		}
		// apply boundaries
		// west
		if (!pinfo->hasNbr(Side<2>::west)) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				double y = pinfo->starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * pinfo->ns[0]] -= efun(pinfo->starts[0], y) * 2 / (h_x * h_x);
			}
		}
		// east
		if (!pinfo->hasNbr(Side<2>::east)) {
			for (int yi = 0; yi < pinfo->ns[1]; yi++) {
				double y = pinfo->starts[1] + h_y / 2.0 + h_y * yi;
				f_vals[yi * pinfo->ns[0] + pinfo->ns[0] - 1]
				-= efun(pinfo->starts[0] + h_x * pinfo->ns[0], y) * 2 / (h_x * h_x);
			}
		}
		// south
		if (!pinfo->hasNbr(Side<2>::south)) {
			for (int xi = 0; xi < pinfo->ns[0]; xi++) {
				double x = pinfo->starts[0] + h_x / 2.0 + h_x * xi;
				f_vals[xi] -= efun(x, pinfo->starts[1]) * 2 / (h_y * h_y);
			}
		}
		// north
		if (!pinfo->hasNbr(Side<2>::north)) {
			for (int xi = 0; xi < pinfo->ns[0]; xi++) {
				double x = pinfo->starts[0] + h_x / 2.0 + h_x * xi;
				f_vals[pinfo->ns[0] * (pinfo->ns[1] - 1) + xi]
				-= efun(x, pinfo->starts[1] + h_y * pinfo->ns[1]) * 2 / (h_y * h_y);
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::fillSolution2d(Domain<2> &domain, Vec u, function<double(double, double, double)> fun,
                          double time)
{
	double *vec;
	VecGetArray(u, &vec);
	for (auto &pinfo : domain.getPatchInfoVector()) {
		double *f = vec + pinfo->local_index * pinfo->ns[0] * pinfo->ns[1];

		double h_x = pinfo->spacings[0];
		double h_y = pinfo->spacings[1];
		for (int yi = 0; yi < pinfo->ns[1]; yi++) {
			for (int xi = 0; xi < pinfo->ns[0]; xi++) {
				double x                  = pinfo->starts[0] + h_x / 2.0 + h_x * xi;
				double y                  = pinfo->starts[1] + h_y / 2.0 + h_y * yi;
				f[yi * pinfo->ns[0] + xi] = fun(x, y, time);
			}
		}
	}
	VecRestoreArray(u, &vec);
}
