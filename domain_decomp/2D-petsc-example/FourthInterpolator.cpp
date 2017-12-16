#include "FourthInterpolator.h"
class Slice
{
	private:
	const double *start;
	int           stride;

	public:
	Slice() {}
	Slice(const double *start, int stride)
	{
		this->start  = start;
		this->stride = stride;
	}

	const double &operator[](size_t idx) { return start[stride * idx]; }
};
Slice getSlice(Domain &d, int n, const double *u_view, Side s)
{
	int   start = d.id_local * n * n;
	Slice retval;
	switch (s) {
		case Side::north:
			retval = Slice(&u_view[start + n * (n - 1)], 1);
			break;
		case Side::east:
			retval = Slice(&u_view[start + (n - 1)], n);
			break;
		case Side::south:
			retval = Slice(&u_view[start], 1);
			break;
		case Side::west:
			retval = Slice(&u_view[start], n);
			break;
	}
	return retval;
}
Slice getInnerSlice(Domain &d, int n, const double *u_view, Side s)
{
	int   start = d.id_local * n * n;
	Slice retval;
	switch (s) {
		case Side::north:
			retval = Slice(&u_view[start + n * (n - 2)], 1);
			break;
		case Side::east:
			retval = Slice(&u_view[start + (n - 2)], n);
			break;
		case Side::south:
			retval = Slice(&u_view[start + n], 1);
			break;
		case Side::west:
			retval = Slice(&u_view[start + 1], n);
			break;
	}
	return retval;
}
void FourthInterpolator::interpolate(Domain &d, const Vec u, Vec interp)
{
	Side s = Side::north;
	do {
		if (d.hasFineNbr(s)) {
			// three sets gammas to take care of
			interpolate(d, s, InterpCase::coarse_from_coarse, u, interp);
			interpolate(d, s, InterpCase::fine_from_coarse_to_fine_on_left, u, interp);
			interpolate(d, s, InterpCase::fine_from_coarse_to_fine_on_right, u, interp);
		} else if (d.hasCoarseNbr(s)) {
			if (d.leftOfCoarse(s)) {
				interpolate(d, s, InterpCase::fine_from_fine_on_left, u, interp);
				interpolate(d, s, InterpCase::coarse_from_fine_on_left, u, interp);
			} else {
				interpolate(d, s, InterpCase::fine_from_fine_on_right, u, interp);
				interpolate(d, s, InterpCase::coarse_from_fine_on_right, u, interp);
			}
		} else if (d.hasNbr(s)) {
			interpolate(d, s, InterpCase::normal, u, interp);
		}
		s++;
	} while (s != Side::north);
}
void FourthInterpolator::interpolate(Domain &d, Side s, InterpCase icase, const Vec u, Vec interp)
{
	int     n = d.n;
	double *interp_view;
	VecGetArray(interp, &interp_view);
	double *u_view;
	VecGetArray(u, &u_view);
	int right  = -1;
	int ctfidx = -1;
	switch (icase) {
		case InterpCase::normal: {
			Slice sl  = getSlice(d, n, u_view, s);
			Slice isl = getInnerSlice(d, n, u_view, s);
			int   idx = n * d.index(s);
			for (int i = 0; i < n; i++) {
				interp_view[idx + i] += 9.0 / 16.0 * sl[i] - 1.0 / 16.0 * isl[i];
			}
		} break;
		case InterpCase::coarse_from_coarse: {
			Slice sl  = getSlice(d, n, u_view, s);
			int   idx = n * d.index(s);
			// beginning case
			interp_view[idx] += ((-sl[2] + 2 * sl[1] - 3 * sl[0]) / 30.0 + sl[0]) / 2;
			// middle cases
			for (int i = 1; i < n - 1; i++) {
				interp_view[idx + i] += (-1.0 / 30 * sl[i - 1] + sl[i] - 1.0 / 30 * sl[i + 1]) / 2;
			}
			// end case
			interp_view[idx + n - 1]
			+= ((-sl[n - 3] + 2 * sl[n - 2] - 3 * sl[n - 1]) / 30.0 + sl[n - 1]) / 2;
		} break;
		case InterpCase::coarse_from_fine_on_left:
			right++;
		case InterpCase::coarse_from_fine_on_right:
			right++;
			{
				Slice sl                = getSlice(d, n, u_view, s);
				Slice isl               = getInnerSlice(d, n, u_view, s);
				int   idx               = n * d.indexCenter(s);
				bool  edge_on_beginning = ((s == Side::north || s == Side::west) && !right)
				                         || ((s == Side::east || s == Side::south) && right);
				if (edge_on_beginning) {
					// middle cases
					for (int i = 0; i < n; i++) {
						interp_view[idx + i / 2] += 1.0 / 6 * sl[i] + 1.0 / 10 * isl[i];
					}
				} else {
					// middle cases
					for (int i = 0; i < n; i++) {
						interp_view[idx + (n + i) / 2] += 1.0 / 6 * sl[i] + 1.0 / 10 * isl[i];
					}
				}
			}
			break;
		case InterpCase::fine_from_fine_on_left:
		case InterpCase::fine_from_fine_on_right: {
			Slice sl  = getSlice(d, n, u_view, s);
			Slice isl = getInnerSlice(d, n, u_view, s);
			int   idx = n * d.index(s);
			// middle cases
			for (int i = 0; i < n; i++) {
				interp_view[idx + i] += 5.0 / 6 * sl[i] - 1.0 / 10 * isl[i];
			}
		} break;
		case InterpCase::fine_from_coarse_to_fine_on_left:
			right++;
			ctfidx = n * d.indexRefinedLeft(s);
		case InterpCase::fine_from_coarse_to_fine_on_right:
			right++;
			if (ctfidx == -1) ctfidx = n * d.indexRefinedRight(s);
			{
				int   idx               = ctfidx;
				Slice sl                = getSlice(d, n, u_view, s);
				bool  edge_on_beginning = ((s == Side::north || s == Side::west) && right)
				                         || ((s == Side::east || s == Side::south) && !right);
				if (edge_on_beginning) {
					// edge case
					interp_view[idx] += (5 * sl[2] - 18 * sl[1] + 45 * sl[0]) / 120.0;
					interp_view[idx + 1] += (-3 * sl[2] + 14 * sl[1] + 21 * sl[0]) / 120.0;
					// middle cases
					for (int i = 2; i < n; i++) {
						if (i % 2 == 0) {
							interp_view[idx + i] += 1.0 / 24 * sl[i / 2 - 1] + 1.0 / 4 * sl[i / 2]
							                        - 1.0 / 40 * sl[i / 2 + 1];
						} else {
							interp_view[idx + i] += 1.0 / 24 * sl[i / 2 + 1] + 1.0 / 4 * sl[i / 2]
							                        - 1.0 / 40 * sl[i / 2 - 1];
						}
					}
				} else {
					// edge case
					interp_view[idx + n - 1]
					+= (5 * sl[n - 3] - 18 * sl[n - 2] + 45 * sl[n - 1]) / 120.0;
					interp_view[idx + n - 2]
					+= (-3 * sl[n - 3] + 14 * sl[n - 2] + 21 * sl[n - 1]) / 120.0;
					// middle cases
					for (int i = n - 3; i >= 0; i--) {
						if ((n + i) % 2 == 0) {
							interp_view[idx + i] += 1.0 / 24 * sl[(n + i) / 2 - 1]
							                        + 1.0 / 4 * sl[(n + i) / 2]
							                        - 1.0 / 40 * sl[(n + i) / 2 + 1];
						} else {
							interp_view[idx + i] += 1.0 / 24 * sl[(n + i) / 2 + 1]
							                        + 1.0 / 4 * sl[(n + i) / 2]
							                        - 1.0 / 40 * sl[(n + i) / 2 - 1];
						}
					}
				}
			}
			break;
	}
	VecRestoreArray(interp, &interp_view);
	VecRestoreArray(u, &u_view);
}
