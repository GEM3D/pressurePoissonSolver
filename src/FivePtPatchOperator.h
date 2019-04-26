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
#include <PatchOperator.h>
class FivePtPatchOperator : public PatchOperator<2>
{
	public:
	void apply(SchurDomain<2> &d, std::shared_ptr<const Vector<2>> u,
	           std::shared_ptr<const Vector<1>> gamma, std::shared_ptr<Vector<2>> f)
	{
		int nx = d.ns[0];
		int ny = d.ns[1];

		double h_x = d.spacings[0];
		double h_y = d.spacings[1];

		LocalData<2>       f_data = f->getLocalData(d.id_local);
		const LocalData<2> u_data = u->getLocalData(d.id_local);

		double center, north, east, south, west;

		// derive in x direction
		// west
		if (d.hasNbr(Side<2>::west)) {
			const LocalData<1> boundary_data
			= gamma->getLocalData(d.getIfaceLocalIndex(Side<2>::west));
			for (int yi = 0; yi < ny; yi++) {
				west            = boundary_data[{yi}];
				center          = u_data[{0, yi}];
				east            = u_data[{1, yi}];
				f_data[{0, yi}] = (2 * west - 3 * center + east) / (h_x * h_x);
			}
		} else if (d.isNeumann(Side<2>::west)) {
			for (int yi = 0; yi < ny; yi++) {
				center          = u_data[{0, yi}];
				east            = u_data[{1, yi}];
				f_data[{0, yi}] = (-center + east) / (h_x * h_x);
			}
		} else {
			for (int yi = 0; yi < ny; yi++) {
				center          = u_data[{0, yi}];
				east            = u_data[{1, yi}];
				f_data[{0, yi}] = (-3 * center + east) / (h_x * h_x);
			}
		}
		// middle
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 1; xi < nx - 1; xi++) {
				west   = u_data[{xi - 1, yi}];
				center = u_data[{xi, yi}];
				east   = u_data[{xi + 1, yi}];

				f_data[{xi, yi}] = (west - 2 * center + east) / (h_x * h_x);
			}
		}
		// east
		if (d.hasNbr(Side<2>::east)) {
			const LocalData<1> boundary_data
			= gamma->getLocalData(d.getIfaceLocalIndex(Side<2>::east));
			for (int yi = 0; yi < ny; yi++) {
				west                 = u_data[{nx - 2, yi}];
				center               = u_data[{nx - 1, yi}];
				east                 = boundary_data[{yi}];
				f_data[{nx - 1, yi}] = (west - 3 * center + 2 * east) / (h_x * h_x);
			}
		} else if (d.isNeumann(Side<2>::east)) {
			for (int yi = 0; yi < ny; yi++) {
				west                 = u_data[{nx - 2, yi}];
				center               = u_data[{nx - 1, yi}];
				f_data[{nx - 1, yi}] = (west - center) / (h_x * h_x);
			}
		} else {
			for (int yi = 0; yi < ny; yi++) {
				west                 = u_data[{nx - 2, yi}];
				center               = u_data[{nx - 1, yi}];
				f_data[{nx - 1, yi}] = (west - 3 * center) / (h_x * h_x);
			}
		}

		// derive in y direction
		// south
		if (d.hasNbr(Side<2>::south)) {
			const LocalData<1> boundary_data
			= gamma->getLocalData(d.getIfaceLocalIndex(Side<2>::south));
			for (int xi = 0; xi < nx; xi++) {
				south  = boundary_data[{xi}];
				center = u_data[{xi, 0}];
				north  = u_data[{xi, 1}];
				f_data[{xi, 0}] += (2 * south - 3 * center + north) / (h_y * h_y);
			}
		} else if (d.isNeumann(Side<2>::south)) {
			for (int xi = 0; xi < nx; xi++) {
				center = u_data[{xi, 0}];
				north  = u_data[{xi, 1}];
				f_data[{xi, 0}] += (-center + north) / (h_y * h_y);
			}
		} else {
			for (int xi = 0; xi < nx; xi++) {
				center = u_data[{xi, 0}];
				north  = u_data[{xi, 1}];
				f_data[{xi, 0}] += (-3 * center + north) / (h_y * h_y);
			}
		}
		// middle
		for (int yi = 1; yi < ny - 1; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u_data[{xi, yi - 1}];
				center = u_data[{xi, yi}];
				north  = u_data[{xi, yi + 1}];

				f_data[{xi, yi}] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
		// north
		if (d.hasNbr(Side<2>::north)) {
			const LocalData<1> boundary_data
			= gamma->getLocalData(d.getIfaceLocalIndex(Side<2>::north));
			for (int xi = 0; xi < nx; xi++) {
				south  = u_data[{xi, ny - 2}];
				center = u_data[{xi, ny - 1}];
				north  = boundary_data[{xi}];
				f_data[{xi, ny - 1}] += (south - 3 * center + 2 * north) / (h_y * h_y);
			}
		} else if (d.isNeumann(Side<2>::north)) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u_data[{xi, ny - 2}];
				center = u_data[{xi, ny - 1}];
				f_data[{xi, ny - 1}] += (south - center) / (h_y * h_y);
			}
		} else {
			for (int xi = 0; xi < nx; xi++) {
				south  = u_data[{xi, nx - 2}];
				center = u_data[{xi, nx - 1}];
				f_data[{xi, nx - 1}] += (south - 3 * center) / (h_y * h_y);
			}
		}
	}
};
#endif
