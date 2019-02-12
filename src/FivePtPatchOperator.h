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

#ifndef FIVEPTPATCHOPERATOR_H
#define FIVEPTPATCHOPERATOR_H
#include "PatchOperator.h"
class FivePtPatchOperator : public PatchOperator<2>
{
	public:
	void apply(SchurDomain<2> &d, const Vec u, const Vec gamma, Vec f)
	{
		int     n     = d.n;
		double  h_x   = d.domain.lengths[0] / n;
		double  h_y   = d.domain.lengths[1] / n;
		int     start = n * n * d.domain.id_local;
		double *u_view, *f_view, *gamma_view;
		VecGetArray(u, &u_view);
		VecGetArray(f, &f_view);
		double *f_ptr = f_view + start;
		double *u_ptr = u_view + start;
		VecGetArray(gamma, &gamma_view);
		const double *boundary_north = nullptr;
		if (d.hasNbr(Side<2>::north)) {
			boundary_north = &gamma_view[d.n * d.getIfaceLocalIndex(Side<2>::north)];
		}
		const double *boundary_east = nullptr;
		if (d.hasNbr(Side<2>::east)) {
			boundary_east = &gamma_view[d.n * d.getIfaceLocalIndex(Side<2>::east)];
		}
		const double *boundary_south = nullptr;
		if (d.hasNbr(Side<2>::south)) {
			boundary_south = &gamma_view[d.n * d.getIfaceLocalIndex(Side<2>::south)];
		}
		const double *boundary_west = nullptr;
		if (d.hasNbr(Side<2>::west)) {
			boundary_west = &gamma_view[d.n * d.getIfaceLocalIndex(Side<2>::west)];
		}
		// integrate in x secton
		double center, north, east, south, west;
		// west
		for (int j = 0; j < n; j++) {
			west = 0;
			if (boundary_west != nullptr) { west = boundary_west[j]; }
			center = u_ptr[j * n];
			east   = u_ptr[j * n + 1];
			if (d.isNeumann(Side<2>::west) && !d.hasNbr(Side<2>::west)) {
				f_ptr[j * n] = (-h_x * west - center + east) / (h_x * h_x);
			} else {
				f_ptr[j * n] = (2 * west - 3 * center + east) / (h_x * h_x);
			}
		}
		// middle
		for (int i = 1; i < n - 1; i++) {
			for (int j = 0; j < n; j++) {
				east   = u_ptr[j * n + i - 1];
				center = u_ptr[j * n + i];
				west   = u_ptr[j * n + i + 1];

				f_ptr[j * n + i] = (west - 2 * center + east) / (h_x * h_x);
			}
		}
		// east
		for (int j = 0; j < n; j++) {
			west   = u_ptr[j * n + n - 2];
			center = u_ptr[j * n + n - 1];
			east   = 0;
			if (boundary_east != nullptr) { east = boundary_east[j]; }
			if (d.isNeumann(Side<2>::east) && !d.hasNbr(Side<2>::east)) {
				f_ptr[j * n + n - 1] = (west - center + h_x * east) / (h_x * h_x);
			} else {
				f_ptr[j * n + n - 1] = (west - 3 * center + 2 * east) / (h_x * h_x);
			}
		}
		// south
		for (int i = 0; i < n; i++) {
			south = 0;
			if (boundary_south != nullptr) { south = boundary_south[i]; }
			center = u_ptr[i];
			north  = u_ptr[n + i];
			if (d.isNeumann(Side<2>::south) && !d.hasNbr(Side<2>::south)) {
				f_ptr[i] += (-h_y * south - center + north) / (h_y * h_y);
			} else {
				f_ptr[i] += (2 * south - 3 * center + north) / (h_y * h_y);
			}
		}
		// middle
		for (int i = 0; i < n; i++) {
			for (int j = 1; j < n - 1; j++) {
				south  = u_ptr[(j - 1) * n + i];
				center = u_ptr[j * n + i];
				north  = u_ptr[(j + 1) * n + i];

				f_ptr[j * n + i] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
		// north
		for (int i = 0; i < n; i++) {
			south  = u_ptr[(n - 2) * n + i];
			center = u_ptr[(n - 1) * n + i];
			north  = 0;
			if (boundary_north != nullptr) { north = boundary_north[i]; }
			if (d.isNeumann(Side<2>::north) && !d.hasNbr(Side<2>::north)) {
				f_ptr[(n - 1) * n + i] += (south - center + h_y * north) / (h_y * h_y);
			} else {
				f_ptr[(n - 1) * n + i] += (south - 3 * center + 2 * north) / (h_y * h_y);
			}
		}
		VecRestoreArray(gamma, &gamma_view);
		VecRestoreArray(u, &u_view);
		VecRestoreArray(f, &f_view);
	}
};
#endif
