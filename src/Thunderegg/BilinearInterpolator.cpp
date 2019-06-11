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
using namespace std;
void BilinearInterpolator::interpolate(const std::vector<SchurDomain<2>> &patches,
                                       std::shared_ptr<const Vector<2>>   u,
                                       std::shared_ptr<Vector<1>>         interp)
{
	for (SchurDomain<2> p : patches) {
		for (Side<2> s : Side<2>::getValues()) {
			if (p.hasNbr(s)) {
				std::deque<int>       idx;
				std::deque<IfaceType> types;

				p.getIfaceInfoPtr(s)->getLocalIndexes(idx);
				p.getIfaceInfoPtr(s)->getIfaceTypes(types);

				for (size_t i = 0; i < idx.size(); i++) {
					interpolate(p, s, idx[i], types[i], u, interp);
				}
			}
		}
	}
}
void BilinearInterpolator::interpolate(SchurDomain<2> &d, std::shared_ptr<const Vector<2>> u,
                                       std::shared_ptr<Vector<1>> interp)
{
	for (Side<2> s : Side<2>::getValues()) {
		if (d.hasNbr(s)) {
			std::deque<int>       idx;
			std::deque<IfaceType> types;

			d.getIfaceInfoPtr(s)->getLocalIndexes(idx);
			d.getIfaceInfoPtr(s)->getIfaceTypes(types);

			for (size_t i = 0; i < idx.size(); i++) {
				interpolate(d, s, idx[i], types[i], u, interp);
			}
		}
	}
}
void BilinearInterpolator::interpolate(SchurDomain<2> &d, Side<2> s, int local_index,
                                       IfaceType itype, std::shared_ptr<const Vector<2>> u,
                                       std::shared_ptr<Vector<1>> interp)
{
	int n = d.ns[!s.axis()];

	LocalData<1>       interp_data = interp->getLocalData(local_index);
	const LocalData<1> sl          = u->getLocalData(d.local_index).getSliceOnSide(s);

	switch (itype.toInt()) {
		case IfaceType::normal: {
			for (int i = 0; i < n; i++) {
				interp_data[{i}] += 0.5 * sl[{i}];
			}
		} break;
		case IfaceType::coarse_to_coarse: {
			// middle cases
			for (int i = 0; i < n; i++) {
				interp_data[{i}] += 1.0 / 3 * sl[{i}];
			}
		} break;
		case IfaceType::fine_to_coarse: {
			if (itype.getOrthant() == 0) {
				// middle cases
				for (int i = 0; i < n; i += 2) {
					interp_data[{i / 2}] += 1.0 / 3 * sl[{i}] + 1.0 / 3 * sl[{i + 1}];
				}
			} else {
				// middle cases
				for (int i = 0; i < n; i += 2) {
					interp_data[{(n + i) / 2}] += 1.0 / 3 * sl[{i}] + 1.0 / 3 * sl[{i + 1}];
				}
			}
		} break;
		case IfaceType::fine_to_fine: {
			// middle cases
			for (int i = 0; i < n; i += 2) {
				interp_data[{i}] += 5.0 / 6 * sl[{i}] - 1.0 / 6 * sl[{i + 1}];
			}
			for (int i = 1; i < n; i += 2) {
				interp_data[{i}] += 5.0 / 6 * sl[{i}] - 1.0 / 6 * sl[{i - 1}];
			}
		} break;
		case IfaceType::coarse_to_fine: {
			if (itype.getOrthant() == 0) {
				for (int i = 0; i < n; i++) {
					interp_data[{i}] += 2.0 / 6 * sl[{i / 2}];
				}
			} else {
				// middle cases
				for (int i = 0; i < n; i++) {
					interp_data[{i}] += 2.0 / 6 * sl[{(n + i) / 2}];
				}
			}
		} break;
	}
}
