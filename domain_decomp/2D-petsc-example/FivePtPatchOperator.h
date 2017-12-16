#ifndef FIVEPTPATCHOPERATOR_H
#define FIVEPTPATCHOPERATOR_H
#include "PatchOperator.h"
class FivePtPatchOperator : public PatchOperator
{
	public:
	void apply(Domain &d, const vector_type &u, const vector_type &gamma, vector_type &f)
	{
		int           n              = d.n;
		double        h_x            = d.x_length / n;
		double        h_y            = d.y_length / n;
		int           start          = n * n * d.id_local;
		auto          u_ptr          = &(u.get1dView()[start]);
		auto          f_ptr          = &(f.get1dViewNonConst()[start]);
		auto          gamma_view     = gamma.get1dView();
		const double *boundary_north = nullptr;
		if (d.hasNbr(Side::north)) {
			boundary_north = &gamma_view[d.n * d.index(Side::north)];
		}
		const double *boundary_east = nullptr;
		if (d.hasNbr(Side::east)) {
			boundary_east = &gamma_view[d.n * d.index(Side::east)];
		}
		const double *boundary_south = nullptr;
		if (d.hasNbr(Side::south)) {
			boundary_south = &gamma_view[d.n * d.index(Side::south)];
		}
		const double *boundary_west = nullptr;
		if (d.hasNbr(Side::west)) {
			boundary_west = &gamma_view[d.n * d.index(Side::west)];
		}
		// integrate in x secton
		double center, north, east, south, west;
		// west
		for (int j = 0; j < n; j++) {
			west = 0;
			if (boundary_west != nullptr) {
				west = boundary_west[j];
			}
			center = u_ptr[j * n];
			east   = u_ptr[j * n + 1];
			if (d.isNeumann(Side::west) && !d.hasNbr(Side::west)) {
				f_ptr[j * n] = (-h_x * west - center + east) / (h_x * h_x);
			} else {
				f_ptr[j * n] = (2 * west - 3 * center + east) / (h_x * h_x);
			}
		}
		// middle
		for (int i = 1; i < n - 1; i++) {
			for (int j = 0; j < n; j++) {
				east   = u_ptr[j * n + i - 1];
				center = u_ptr[j * n + i];
				west   = u_ptr[j * n + i + 1];

				f_ptr[j * n + i] = (west - 2 * center + east) / (h_x * h_x);
			}
		}
		// east
		for (int j = 0; j < n; j++) {
			west   = u_ptr[j * n + n - 2];
			center = u_ptr[j * n + n - 1];
			east   = 0;
			if (boundary_east != nullptr) {
				east = boundary_east[j];
			}
			if (d.isNeumann(Side::east) && !d.hasNbr(Side::east)) {
				f_ptr[j * n + n - 1] = (west - center + h_x * east) / (h_x * h_x);
			} else {
				f_ptr[j * n + n - 1] = (west - 3 * center + 2 * east) / (h_x * h_x);
			}
		}
		// south
		for (int i = 0; i < n; i++) {
			south = 0;
			if (boundary_south != nullptr) {
				south = boundary_south[i];
			}
			center = u_ptr[i];
			north  = u_ptr[n + i];
			if (d.isNeumann(Side::south) && !d.hasNbr(Side::south)) {
				f_ptr[i] += (-h_y * south - center + north) / (h_y * h_y);
			} else {
				f_ptr[i] += (2 * south - 3 * center + north) / (h_y * h_y);
			}
		}
		// middle
		for (int i = 0; i < n; i++) {
			for (int j = 1; j < n - 1; j++) {
				south  = u_ptr[(j - 1) * n + i];
				center = u_ptr[j * n + i];
				north  = u_ptr[(j + 1) * n + i];

				f_ptr[j * n + i] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
		// north
		for (int i = 0; i < n; i++) {
			south  = u_ptr[(n - 2) * n + i];
			center = u_ptr[(n - 1) * n + i];
			north  = 0;
			if (boundary_north != nullptr) {
				north = boundary_north[i];
			}
			if (d.isNeumann(Side::north) && !d.hasNbr(Side::north)) {
				f_ptr[(n - 1) * n + i] += (south - center + h_y * north) / (h_y * h_y);
			} else {
				f_ptr[(n - 1) * n + i] += (south - 3 * center + 2 * north) / (h_y * h_y);
			}
		}
	}
};
#endif
