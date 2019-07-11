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
#include <Thunderegg/PatchOperator.h>
class FivePtPatchOperator : public PatchOperator<2>
{
	public:
	void applyWithInterface(SchurInfo<2> &sinfo, const LocalData<2> u,
	                        std::shared_ptr<const Vector<1>> gamma, LocalData<2> f) override
	{
		int nx = sinfo.pinfo->ns[0];
		int ny = sinfo.pinfo->ns[1];

		double h_x = sinfo.pinfo->spacings[0];
		double h_y = sinfo.pinfo->spacings[1];

		double center, north, east, south, west;

		// derive in x direction
		// west
		if (sinfo.pinfo->hasNbr(Side<2>::west)) {
			const LocalData<1> boundary_data
			= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<2>::west));
			for (int yi = 0; yi < ny; yi++) {
				west       = boundary_data[{yi}];
				center     = u[{0, yi}];
				east       = u[{1, yi}];
				f[{0, yi}] = (2 * west - 3 * center + east) / (h_x * h_x);
			}
		} else if (sinfo.pinfo->isNeumann(Side<2>::west)) {
			for (int yi = 0; yi < ny; yi++) {
				center     = u[{0, yi}];
				east       = u[{1, yi}];
				f[{0, yi}] = (-center + east) / (h_x * h_x);
			}
		} else {
			for (int yi = 0; yi < ny; yi++) {
				center     = u[{0, yi}];
				east       = u[{1, yi}];
				f[{0, yi}] = (-3 * center + east) / (h_x * h_x);
			}
		}
		// middle
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 1; xi < nx - 1; xi++) {
				west   = u[{xi - 1, yi}];
				center = u[{xi, yi}];
				east   = u[{xi + 1, yi}];

				f[{xi, yi}] = (west - 2 * center + east) / (h_x * h_x);
			}
		}
		// east
		if (sinfo.pinfo->hasNbr(Side<2>::east)) {
			const LocalData<1> boundary_data
			= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<2>::east));
			for (int yi = 0; yi < ny; yi++) {
				west            = u[{nx - 2, yi}];
				center          = u[{nx - 1, yi}];
				east            = boundary_data[{yi}];
				f[{nx - 1, yi}] = (west - 3 * center + 2 * east) / (h_x * h_x);
			}
		} else if (sinfo.pinfo->isNeumann(Side<2>::east)) {
			for (int yi = 0; yi < ny; yi++) {
				west            = u[{nx - 2, yi}];
				center          = u[{nx - 1, yi}];
				f[{nx - 1, yi}] = (west - center) / (h_x * h_x);
			}
		} else {
			for (int yi = 0; yi < ny; yi++) {
				west            = u[{nx - 2, yi}];
				center          = u[{nx - 1, yi}];
				f[{nx - 1, yi}] = (west - 3 * center) / (h_x * h_x);
			}
		}

		// derive in y direction
		// south
		if (sinfo.pinfo->hasNbr(Side<2>::south)) {
			const LocalData<1> boundary_data
			= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<2>::south));
			for (int xi = 0; xi < nx; xi++) {
				south  = boundary_data[{xi}];
				center = u[{xi, 0}];
				north  = u[{xi, 1}];
				f[{xi, 0}] += (2 * south - 3 * center + north) / (h_y * h_y);
			}
		} else if (sinfo.pinfo->isNeumann(Side<2>::south)) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, 0}];
				north  = u[{xi, 1}];
				f[{xi, 0}] += (-center + north) / (h_y * h_y);
			}
		} else {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, 0}];
				north  = u[{xi, 1}];
				f[{xi, 0}] += (-3 * center + north) / (h_y * h_y);
			}
		}
		// middle
		for (int yi = 1; yi < ny - 1; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, yi - 1}];
				center = u[{xi, yi}];
				north  = u[{xi, yi + 1}];

				f[{xi, yi}] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
		// north
		if (sinfo.pinfo->hasNbr(Side<2>::north)) {
			const LocalData<1> boundary_data
			= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<2>::north));
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, ny - 2}];
				center = u[{xi, ny - 1}];
				north  = boundary_data[{xi}];
				f[{xi, ny - 1}] += (south - 3 * center + 2 * north) / (h_y * h_y);
			}
		} else if (sinfo.pinfo->isNeumann(Side<2>::north)) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, ny - 2}];
				center = u[{xi, ny - 1}];
				f[{xi, ny - 1}] += (south - center) / (h_y * h_y);
			}
		} else {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, nx - 2}];
				center = u[{xi, nx - 1}];
				f[{xi, nx - 1}] += (south - 3 * center) / (h_y * h_y);
			}
		}
	}
	void addInterfaceToRHS(SchurInfo<2> &sinfo, std::shared_ptr<const Vector<1>> gamma,
	                       LocalData<2> f) override
	{
		for (Side<2> s : Side<2>::getValues()) {
			if (sinfo.pinfo->hasNbr(s)) {
				const LocalData<1> gamma_view = gamma->getLocalData(sinfo.getIfaceLocalIndex(s));

				LocalData<1> slice = f.getSliceOnSide(s);

				double h2 = pow(sinfo.pinfo->spacings[s.axis()], 2);

				nested_loop<1>(
				gamma_view.getStart(), gamma_view.getEnd(),
				[&](std::array<int, 1> coord) { slice[coord] -= 2.0 / h2 * gamma_view[coord]; });
			}
		}
	}
	void apply(SchurInfo<2> &sinfo, const LocalData<2> u, LocalData<2> f) override
	{
		int nx = sinfo.pinfo->ns[0];
		int ny = sinfo.pinfo->ns[1];

		double h_x = sinfo.pinfo->spacings[0];
		double h_y = sinfo.pinfo->spacings[1];

		double center, north, east, south, west;

		// derive in x direction
		// west
		if (sinfo.pinfo->isNeumann(Side<2>::west)) {
			for (int yi = 0; yi < ny; yi++) {
				center     = u[{0, yi}];
				east       = u[{1, yi}];
				f[{0, yi}] = (-center + east) / (h_x * h_x);
			}
		} else {
			for (int yi = 0; yi < ny; yi++) {
				center     = u[{0, yi}];
				east       = u[{1, yi}];
				f[{0, yi}] = (-3 * center + east) / (h_x * h_x);
			}
		}
		// middle
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 1; xi < nx - 1; xi++) {
				west   = u[{xi - 1, yi}];
				center = u[{xi, yi}];
				east   = u[{xi + 1, yi}];

				f[{xi, yi}] = (west - 2 * center + east) / (h_x * h_x);
			}
		}
		// east
		if (sinfo.pinfo->isNeumann(Side<2>::east)) {
			for (int yi = 0; yi < ny; yi++) {
				west            = u[{nx - 2, yi}];
				center          = u[{nx - 1, yi}];
				f[{nx - 1, yi}] = (west - center) / (h_x * h_x);
			}
		} else {
			for (int yi = 0; yi < ny; yi++) {
				west            = u[{nx - 2, yi}];
				center          = u[{nx - 1, yi}];
				f[{nx - 1, yi}] = (west - 3 * center) / (h_x * h_x);
			}
		}

		// derive in y direction
		// south
		if (sinfo.pinfo->isNeumann(Side<2>::south)) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, 0}];
				north  = u[{xi, 1}];
				f[{xi, 0}] += (-center + north) / (h_y * h_y);
			}
		} else {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, 0}];
				north  = u[{xi, 1}];
				f[{xi, 0}] += (-3 * center + north) / (h_y * h_y);
			}
		}
		// middle
		for (int yi = 1; yi < ny - 1; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, yi - 1}];
				center = u[{xi, yi}];
				north  = u[{xi, yi + 1}];

				f[{xi, yi}] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
		// north
		if (sinfo.pinfo->isNeumann(Side<2>::north)) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, ny - 2}];
				center = u[{xi, ny - 1}];
				f[{xi, ny - 1}] += (south - center) / (h_y * h_y);
			}
		} else {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, nx - 2}];
				center = u[{xi, nx - 1}];
				f[{xi, nx - 1}] += (south - 3 * center) / (h_y * h_y);
			}
		}
	}
};
#endif
