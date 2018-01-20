#include "TriLinInterp.h"
#include "Utils.h"
using namespace Utils;
void TriLinInterp::interpolate(SchurDomain &d, const Vec u, Vec interp)
{
	for (Side s : getSideValues()) {
		if (d.hasNbr(s)) {
			interpolate(d, s, d.getIfaceLocalIndex(s), IfaceType::normal, u, interp);
		}
	}
}
void TriLinInterp::interpolate(SchurDomain &d, Side s, int local_index, IfaceType itype,
                               const Vec u, Vec interp)
{
	int     n = d.n;
	double *interp_view;
	VecGetArray(interp, &interp_view);
	double *u_view;
	VecGetArray(u, &u_view);
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
		default:
			break;
	}
	VecRestoreArray(interp, &interp_view);
	VecRestoreArray(u, &u_view);
}
