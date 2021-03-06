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

#include "TriLinInterp.h"
void TriLinInterp::interpolate(const std::vector<SchurInfo<3>> &patches,
                               std::shared_ptr<const Vector<3>> u,
                               std::shared_ptr<Vector<2>>       interp)
{
	for (SchurInfo<3> p : patches) {
		for (Side<3> s : Side<3>::getValues()) {
			if (p.pinfo->hasNbr(s)) {
				std::deque<int>          idx;
				std::deque<IfaceType<3>> types;

				p.getIfaceInfoPtr(s)->getLocalIndexes(idx);
				p.getIfaceInfoPtr(s)->getIfaceTypes(types);

				for (size_t i = 0; i < idx.size(); i++) {
					interpolate(p, s, idx[i], types[i], u, interp);
				}
			}
		}
	}
}
void TriLinInterp::interpolate(SchurInfo<3> &sinfo, std::shared_ptr<const Vector<3>> u,
                               std::shared_ptr<Vector<2>> interp)
{
	for (Side<3> s : Side<3>::getValues()) {
		if (sinfo.pinfo->hasNbr(s)) {
			std::deque<int>          idx;
			std::deque<IfaceType<3>> types;

			sinfo.getIfaceInfoPtr(s)->getLocalIndexes(idx);
			sinfo.getIfaceInfoPtr(s)->getIfaceTypes(types);

			for (size_t i = 0; i < idx.size(); i++) {
				interpolate(sinfo, s, idx[i], types[i], u, interp);
			}
		}
	}
}
void TriLinInterp::interpolate(SchurInfo<3> &sinfo, Side<3> s, int local_index, IfaceType<3> itype,
                               std::shared_ptr<const Vector<3>> u,
                               std::shared_ptr<Vector<2>>       interp)
{
	std::array<int, 2> ns;
	for (int i = 0; i < s.axis(); i++) {
		ns[i] = sinfo.pinfo->ns[i];
	}
	for (int i = s.axis(); i < 2; i++) {
		ns[i] = sinfo.pinfo->ns[i + 1];
	}
	int nx = ns[0];
	int ny = ns[1];

	LocalData<2>       interp_data = interp->getLocalData(local_index);
	const LocalData<2> sl          = u->getLocalData(sinfo.pinfo->local_index).getSliceOnSide(s);

	switch (itype.toInt()) {
		case IfaceType<3>::normal: {
			for (int yi = 0; yi < ny; yi++) {
				for (int xi = 0; xi < nx; xi++) {
					interp_data[{xi, yi}] += 0.5 * sl[{xi, yi}];
				}
			}
		} break;
		case IfaceType<3>::fine_to_fine: {
			for (int yi = 0; yi < ny / 2; yi++) {
				for (int xi = 0; xi < nx / 2; xi++) {
					double a     = sl[{xi * 2, yi * 2}];
					double b     = sl[{xi * 2 + 1, yi * 2}];
					double c     = sl[{xi * 2, yi * 2 + 1}];
					double sinfo = sl[{xi * 2 + 1, yi * 2 + 1}];
					interp_data[{xi * 2, yi * 2}] += (11 * a - b - c - sinfo) / 12.0;
					interp_data[{xi * 2 + 1, yi * 2}] += (-a + 11 * b - c - sinfo) / 12.0;
					interp_data[{xi * 2, yi * 2 + 1}] += (-a - b + 11 * c - sinfo) / 12.0;
					interp_data[{xi * 2 + 1, yi * 2 + 1}] += (-a - b - c + 11 * sinfo) / 12.0;
				}
			}
		} break;
		case IfaceType<3>::coarse_to_fine:
			switch (itype.getOrthant().toInt()) {
				case 0: {
					for (int yi = 0; yi < ny; yi++) {
						for (int xi = 0; xi < nx; xi++) {
							interp_data[{xi, yi}] += 4.0 * sl[{(xi) / 2, (yi) / 2}] / 12.0;
						}
					}
				} break;
				case 1: {
					for (int yi = 0; yi < ny; yi++) {
						for (int xi = 0; xi < nx; xi++) {
							interp_data[{xi, yi}] += 4.0 * sl[{(xi + nx) / 2, (yi) / 2}] / 12.0;
						}
					}
				} break;
				case 2: {
					for (int yi = 0; yi < ny; yi++) {
						for (int xi = 0; xi < nx; xi++) {
							interp_data[{xi, yi}] += 4.0 * sl[{(xi) / 2, (yi + ny) / 2}] / 12.0;
						}
					}
				} break;
				case 3: {
					for (int yi = 0; yi < ny; yi++) {
						for (int xi = 0; xi < nx; xi++) {
							interp_data[{xi, yi}]
							+= 4.0 * sl[{(xi + nx) / 2, (yi + ny) / 2}] / 12.0;
						}
					}
				} break;
			}
			break;
		case IfaceType<3>::coarse_to_coarse: {
			for (int yi = 0; yi < ny; yi++) {
				for (int xi = 0; xi < nx; xi++) {
					interp_data[{xi, yi}] += 2.0 / 6.0 * sl[{xi, yi}];
				}
			}
		} break;
		case IfaceType<3>::fine_to_coarse:
			switch (itype.getOrthant().toInt()) {
				case 0: {
					for (int yi = 0; yi < ny; yi++) {
						for (int xi = 0; xi < nx; xi++) {
							interp_data[{(xi) / 2, (yi) / 2}] += 1.0 / 6.0 * sl[{xi, yi}];
						}
					}
				} break;
				case 1: {
					for (int yi = 0; yi < ny; yi++) {
						for (int xi = 0; xi < nx; xi++) {
							interp_data[{(xi + nx) / 2, (yi) / 2}] += 1.0 / 6.0 * sl[{xi, yi}];
						}
					}
				} break;
				case 2: {
					for (int yi = 0; yi < ny; yi++) {
						for (int xi = 0; xi < nx; xi++) {
							interp_data[{(xi) / 2, (yi + ny) / 2}] += 1.0 / 6.0 * sl[{xi, yi}];
						}
					}
				} break;
				case 3: {
					for (int yi = 0; yi < ny; yi++) {
						for (int xi = 0; xi < nx; xi++) {
							interp_data[{(xi + nx) / 2, (yi + ny) / 2}] += 1.0 / 6.0 * sl[{xi, yi}];
						}
					}
				} break;
			}
			break;
	}
}
