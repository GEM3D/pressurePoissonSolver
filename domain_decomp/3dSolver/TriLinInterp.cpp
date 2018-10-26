#include "TriLinInterp.h"
#include "Utils.h"
using namespace Utils;
void TriLinInterp::interpolate(SchurDomain<3> &d, const Vec u, Vec interp)
{
	for (Side<3> s : Side<3>::getValues()) {
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
void TriLinInterp::interpolate(SchurDomain<3> &d, Side<3> s, int local_index, IfaceType itype,
                               const Vec u, Vec interp)
{
	int     n = d.n;
	double *interp_view;
	VecGetArray(interp, &interp_view);
	double *      u_view;
	const double *const_u_view;
	VecGetArrayRead(u, &const_u_view);
	u_view  = const_cast<double *>(const_u_view);
	int idx = local_index * n * n;
	switch (itype) {
		case IfaceType::normal: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + xi + yi * n] += 0.5 * sl(xi, yi);
				}
			}
		} break;
		case IfaceType::fine_to_fine_0:
		case IfaceType::fine_to_fine_1:
		case IfaceType::fine_to_fine_2:
		case IfaceType::fine_to_fine_3: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n / 2; yi++) {
				for (int xi = 0; xi < n / 2; xi++) {
					double a = sl(xi * 2, yi * 2);
					double b = sl(xi * 2 + 1, yi * 2);
					double c = sl(xi * 2, yi * 2 + 1);
					double d = sl(xi * 2 + 1, yi * 2 + 1);
					interp_view[idx + (xi * 2) + (yi * 2) * n] += (11 * a - b - c - d) / 12.0;
					interp_view[idx + (xi * 2 + 1) + (yi * 2) * n] += (-a + 11 * b - c - d) / 12.0;
					interp_view[idx + (xi * 2) + (yi * 2 + 1) * n] += (-a - b + 11 * c - d) / 12.0;
					interp_view[idx + (xi * 2 + 1) + (yi * 2 + 1) * n]
					+= (-a - b - c + 11 * d) / 12.0;
				}
			}
		} break;
		case IfaceType::coarse_to_fine_0: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + xi + yi * n] += 4.0 * sl((xi) / 2, (yi) / 2) / 12.0;
				}
			}
		} break;
		case IfaceType::coarse_to_fine_1: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + xi + yi * n] += 4.0 * sl((xi + n) / 2, (yi) / 2) / 12.0;
				}
			}
		} break;
		case IfaceType::coarse_to_fine_2: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + xi + yi * n] += 4.0 * sl((xi) / 2, (yi + n) / 2) / 12.0;
				}
			}
		} break;
		case IfaceType::coarse_to_fine_3: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + xi + yi * n] += 4.0 * sl((xi + n) / 2, (yi + n) / 2) / 12.0;
				}
			}
		} break;
		case IfaceType::coarse_to_coarse: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + xi + yi * n] += 2.0 / 6.0 * sl(xi, yi);
				}
			}
		} break;
		case IfaceType::fine_to_coarse_0: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + (xi) / 2 + (yi) / 2 * n] += 1.0 / 6.0 * sl(xi, yi);
				}
			}
		} break;
		case IfaceType::fine_to_coarse_1: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + (xi + n) / 2 + (yi) / 2 * n] += 1.0 / 6.0 * sl(xi, yi);
				}
			}
		} break;
		case IfaceType::fine_to_coarse_2: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + (xi) / 2 + (yi + n) / 2 * n] += 1.0 / 6.0 * sl(xi, yi);
				}
			}
		} break;
		case IfaceType::fine_to_coarse_3: {
			Slice sl = getSlice(d, u_view, s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + (xi + n) / 2 + (yi + n) / 2 * n] += 1.0 / 6.0 * sl(xi, yi);
				}
			}
		} break;
	}
	VecRestoreArray(interp, &interp_view);
	VecRestoreArrayRead(u, &const_u_view);
}
