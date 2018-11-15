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
	int right  = -1;
	int ctfidx = -1;
    int idx = local_index*n;
	switch (itype.toInt()) {
		case IfaceType::normal: {
			Slice<1> sl  = getSlice<2>(d, u_view, s);
			for (int i = 0; i < n; i++) {
				interp_view[idx + i] += 0.5 * sl({i});
			}
		} break;
                                /*
		case IfaceType::coarse_from_coarse: {
			Slice sl  = getSlice(d, n, u_view, s);
			// middle cases
			for (int i = 0; i < n; i++) {
				interp_view[idx + i] += 1.0 / 3 * sl[i];
			}
		} break;
		case IfaceType::coarse_from_fine_on_left:
			right++;
		case IfaceType::coarse_from_fine_on_right:
			right++;
			{
				Slice sl                = getSlice(d, n, u_view, s);
				bool  edge_on_beginning = ((s == Side::north || s == Side::west) && !right)
				                         || ((s == Side::east || s == Side::south) && right);
				if (edge_on_beginning) {
					// middle cases
					for (int i = 0; i < n; i += 2) {
						interp_view[idx + i / 2] += 1.0 / 3 * sl[i] + 1.0 / 3 * sl[i + 1];
					}
				} else {
					// middle cases
					for (int i = 0; i < n; i += 2) {
						interp_view[idx + (n + i) / 2] += 1.0 / 3 * sl[i] + 1.0 / 3 * sl[i + 1];
					}
				}
			}
			break;
		case IfaceType::fine_from_fine: {
			Slice sl  = getSlice(d, n, u_view, s);
			int   idx = n * d.index(s);
			// middle cases
			for (int i = 0; i < n; i += 2) {
				interp_view[idx + i] += 5.0 / 6 * sl[i] - 1.0 / 6 * sl[i + 1];
			}
			for (int i = 1; i < n; i += 2) {
				interp_view[idx + i] += 5.0 / 6 * sl[i] - 1.0 / 6 * sl[i - 1];
			}
		} break;
		case IfaceType::fine_from_coarse_to_fine_on_right: {
			Slice sl                = getSlice(d, n, u_view, s);
			bool  edge_on_beginning = ((s == Side::north || s == Side::west) && right)
			                         || ((s == Side::east || s == Side::south) && !right);
			if (edge_on_beginning) {
				for (int i = 0; i < n; i++) {
					interp_view[idx + i] += 2.0 / 6 * sl[i / 2];
				}
			} else {
				// middle cases
				for (int i = 0; i < n; i++) {
					interp_view[idx + i] += 2.0 / 6 * sl[(n + i) / 2];
				}
			}
		} break;
        */
	}

	VecRestoreArray(interp, &interp_view);
	VecRestoreArray(u, &u_view);
}
