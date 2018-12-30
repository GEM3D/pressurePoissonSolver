#ifndef SEVENPTPATCHOPERATOR_H
#define SEVENPTPATCHOPERATOR_H
#include "PatchOperator.h"
#include "Utils.h"
class SevenPtPatchOperator : public PatchOperator<3>
{
	public:
	void apply(SchurDomain<3> &d, const Vec u, const Vec gamma, Vec f)
	{
		using namespace Utils;
		int           n     = d.n;
		double        h_x   = d.domain.lengths[0] / n;
		double        h_y   = d.domain.lengths[1] / n;
		int           start = n * n * n * d.local_index;
		const double *u_view;
		double *      f_view, *gamma_view;
		VecGetArrayRead(u, &u_view);
		VecGetArray(f, &f_view);
		double *      f_ptr = f_view + start;
		const double *u_ptr = u_view + start;
		VecGetArray(gamma, &gamma_view);

		const double *boundary_west = nullptr;
		if (d.hasNbr(Side<3>::west)) {
			boundary_west = &gamma_view[n * n * d.getIfaceLocalIndex(Side<3>::west)];
		}
		const double *boundary_east = nullptr;
		if (d.hasNbr(Side<3>::east)) {
			boundary_east = &gamma_view[n * n * d.getIfaceLocalIndex(Side<3>::east)];
		}
		const double *boundary_south = nullptr;
		if (d.hasNbr(Side<3>::south)) {
			boundary_south = &gamma_view[n * n * d.getIfaceLocalIndex(Side<3>::south)];
		}
		const double *boundary_north = nullptr;
		if (d.hasNbr(Side<3>::north)) {
			boundary_north = &gamma_view[n * n * d.getIfaceLocalIndex(Side<3>::north)];
		}
		const double *boundary_bottom = nullptr;
		if (d.hasNbr(Side<3>::bottom)) {
			boundary_bottom = &gamma_view[n * n * d.getIfaceLocalIndex(Side<3>::bottom)];
		}
		const double *boundary_top = nullptr;
		if (d.hasNbr(Side<3>::top)) {
			boundary_top = &gamma_view[n * n * d.getIfaceLocalIndex(Side<3>::top)];
		}

		double center, north, east, south, west, bottom, top;

		// derive in x direction
		// west
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				west = 0;
				if (boundary_west != nullptr) { west = boundary_west[yi + zi * n]; }
				center = u_ptr[index(n, 0, yi, zi)];
				east   = u_ptr[index(n, 1, yi, zi)];
				if (d.isNeumann(Side<3>::west) && !d.hasNbr(Side<3>::west)) {
					f_ptr[index(n, 0, yi, zi)] = (-h_x * west - center + east) / (h_x * h_x);
				} else {
					f_ptr[index(n, 0, yi, zi)] = (2 * west - 3 * center + east) / (h_x * h_x);
				}
			}
		}
		// middle
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 1; xi < n - 1; xi++) {
					west   = u_ptr[index(n, xi - 1, yi, zi)];
					center = u_ptr[index(n, xi, yi, zi)];
					east   = u_ptr[index(n, xi + 1, yi, zi)];

					f_ptr[index(n, xi, yi, zi)] = (west - 2 * center + east) / (h_x * h_x);
				}
			}
		}
		// east
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				west   = u_ptr[index(n, n - 2, yi, zi)];
				center = u_ptr[index(n, n - 1, yi, zi)];
				east   = 0;
				if (boundary_east != nullptr) { east = boundary_east[yi + zi * n]; }
				if (d.isNeumann(Side<3>::east) && !d.hasNbr(Side<3>::east)) {
					f_ptr[index(n, n - 1, yi, zi)] = (west - center + h_x * east) / (h_x * h_x);
				} else {
					f_ptr[index(n, n - 1, yi, zi)] = (west - 3 * center + 2 * east) / (h_x * h_x);
				}
			}
		}

		// derive in y direction
		// south
		for (int zi = 0; zi < n; zi++) {
			for (int xi = 0; xi < n; xi++) {
				south = 0;
				if (boundary_south != nullptr) { south = boundary_south[xi + zi * n]; }
				center = u_ptr[index(n, xi, 0, zi)];
				north  = u_ptr[index(n, xi, 1, zi)];
				if (d.isNeumann(Side<3>::south) && !d.hasNbr(Side<3>::south)) {
					f_ptr[index(n, xi, 0, zi)] += (-h_y * south - center + north) / (h_y * h_y);
				} else {
					f_ptr[index(n, xi, 0, zi)] += (2 * south - 3 * center + north) / (h_y * h_y);
				}
			}
		}
		// middle
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 1; yi < n - 1; yi++) {
				for (int xi = 0; xi < n; xi++) {
					south  = u_ptr[index(n, xi, yi - 1, zi)];
					center = u_ptr[index(n, xi, yi, zi)];
					north  = u_ptr[index(n, xi, yi + 1, zi)];

					f_ptr[index(n, xi, yi, zi)] += (south - 2 * center + north) / (h_y * h_y);
				}
			}
		}
		// north
		for (int zi = 0; zi < n; zi++) {
			for (int xi = 0; xi < n; xi++) {
				south  = u_ptr[index(n, xi, n - 2, zi)];
				center = u_ptr[index(n, xi, n - 1, zi)];
				north  = 0;
				if (boundary_north != nullptr) { north = boundary_north[xi + zi * n]; }
				if (d.isNeumann(Side<3>::north) && !d.hasNbr(Side<3>::north)) {
					f_ptr[index(n, xi, n - 1, zi)] += (south - center + h_y * north) / (h_y * h_y);
				} else {
					f_ptr[index(n, xi, n - 1, zi)]
					+= (south - 3 * center + 2 * north) / (h_y * h_y);
				}
			}
		}

		// derive in z direction
		// bottom
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				bottom = 0;
				if (boundary_bottom != nullptr) { bottom = boundary_bottom[xi + yi * n]; }
				center = u_ptr[index(n, xi, yi, 0)];
				top    = u_ptr[index(n, xi, yi, 1)];
				if (d.isNeumann(Side<3>::bottom) && !d.hasNbr(Side<3>::bottom)) {
					f_ptr[index(n, xi, yi, 0)] += (-h_y * bottom - center + top) / (h_y * h_y);
				} else {
					f_ptr[index(n, xi, yi, 0)] += (2 * bottom - 3 * center + top) / (h_y * h_y);
				}
			}
		}
		// middle
		for (int zi = 1; zi < n - 1; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					bottom = u_ptr[index(n, xi, yi, zi - 1)];
					center = u_ptr[index(n, xi, yi, zi)];
					top    = u_ptr[index(n, xi, yi, zi + 1)];

					f_ptr[index(n, xi, yi, zi)] += (bottom - 2 * center + top) / (h_y * h_y);
				}
			}
		}
		// top
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				bottom = u_ptr[index(n, xi, yi, n - 2)];
				center = u_ptr[index(n, xi, yi, n - 1)];
				top    = 0;
				if (boundary_top != nullptr) { top = boundary_top[xi + yi * n]; }
				if (d.isNeumann(Side<3>::top) && !d.hasNbr(Side<3>::top)) {
					f_ptr[index(n, xi, yi, n - 1)] += (bottom - center + h_y * top) / (h_y * h_y);
				} else {
					f_ptr[index(n, xi, yi, n - 1)] += (bottom - 3 * center + 2 * top) / (h_y * h_y);
				}
			}
		}
		VecRestoreArray(gamma, &gamma_view);
		VecRestoreArrayRead(u, &u_view);
		VecRestoreArray(f, &f_view);
	}
};
#endif
