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

#include "BilinearInterpolator.h"
#include "Utils.h"
using namespace std;
using namespace Utils;
void BilinearInterpolator::interpolate(SchurDomain<2> &d, const Vec u, Vec interp)
{
	for (Side<2> s : Side<2>::getValues()) {
		if (d.hasNbr(s)) {
			std::deque<int>       idx;
			std::deque<IfaceType> types;
			d.getIfaceInfoPtr(s)->getIdxAndTypes(idx, types);
			for (size_t i = 0; i < idx.size(); i++) {
				interpolate(d, s, idx[i], types[i], u, interp);
			}
		}
	}
}
void BilinearInterpolator::interpolate(SchurDomain<2> &d, Side<2> s, int local_index,
                                       IfaceType itype, const Vec u, Vec interp)
{
	int     n = d.n;
	double *interp_view;
	VecGetArray(interp, &interp_view);
	double *u_view;
	VecGetArray(u, &u_view);
	int idx = local_index * n;
	switch (itype.toInt()) {
		case IfaceType::normal: {
			Slice<1> sl = getSlice<2>(d, u_view, s);
			for (int i = 0; i < n; i++) {
				interp_view[idx + i] += 0.5 * sl({i});
			}
		} break;
		case IfaceType::coarse_to_coarse: {
			Slice<1> sl = getSlice<2>(d, u_view, s);
			// middle cases
			for (int i = 0; i < n; i++) {
				interp_view[idx + i] += 1.0 / 3 * sl({i});
			}
		} break;
		case IfaceType::fine_to_coarse: {
			Slice<1> sl = getSlice<2>(d, u_view, s);
			if (itype.getOrthant() == 0) {
				// middle cases
				for (int i = 0; i < n; i += 2) {
					interp_view[idx + i / 2] += 1.0 / 3 * sl({i}) + 1.0 / 3 * sl({i + 1});
				}
			} else {
				// middle cases
				for (int i = 0; i < n; i += 2) {
					interp_view[idx + (n + i) / 2] += 1.0 / 3 * sl({i}) + 1.0 / 3 * sl({i + 1});
				}
			}
		} break;
		case IfaceType::fine_to_fine: {
			Slice<1> sl = getSlice<2>(d, u_view, s);
			// middle cases
			for (int i = 0; i < n; i += 2) {
				interp_view[idx + i] += 5.0 / 6 * sl({i}) - 1.0 / 6 * sl({i + 1});
			}
			for (int i = 1; i < n; i += 2) {
				interp_view[idx + i] += 5.0 / 6 * sl({i}) - 1.0 / 6 * sl({i - 1});
			}
		} break;
		case IfaceType::coarse_to_fine: {
			Slice<1> sl = getSlice<2>(d, u_view, s);
			if (itype.getOrthant() == 0) {
				for (int i = 0; i < n; i++) {
					interp_view[idx + i] += 2.0 / 6 * sl({i / 2});
				}
			} else {
				// middle cases
				for (int i = 0; i < n; i++) {
					interp_view[idx + i] += 2.0 / 6 * sl({(n + i) / 2});
				}
			}
		} break;
	}

	VecRestoreArray(interp, &interp_view);
	VecRestoreArray(u, &u_view);
}
