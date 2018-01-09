#include "TriLinInterp.h"
#include "Utils.h"
using namespace Utils;
void TriLinInterp::interpolate(Domain &d, const Vec u, Vec interp)
{
	Side s = Side::west;
	do {
		if (d.hasNbr(s)) {
			interpolate(d, s, IfaceType::normal, u, interp);
		}
		s++;
	} while (s != Side::west);
}
void TriLinInterp::interpolate(Domain &d, Side s, IfaceType itype, const Vec u, Vec interp)
{
	int     n = d.n;
	double *interp_view;
	VecGetArray(interp, &interp_view);
	double *u_view;
	VecGetArray(u, &u_view);
	switch (itype) {
		case IfaceType::normal: {
			Slice sl  = getSlice(d, u_view, s);
			int   idx = n * n * d.index(s);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					interp_view[idx + xi + yi * n] += 0.5 * sl(xi, yi);
				}
			}
		} break;
        case IfaceType::coarse_from_coarse:
        case IfaceType::coarse_from_fine_0:
		default:
			break;
	}
	VecRestoreArray(interp, &interp_view);
	VecRestoreArray(u, &u_view);
}
