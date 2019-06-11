#include "SevenPtPatchOperator.h"
void SevenPtPatchOperator::apply(SchurDomain<3> &d, std::shared_ptr<const Vector<3>> u,
                                 std::shared_ptr<const Vector<2>> gamma,
                                 std::shared_ptr<Vector<3>>       f)

{
	int nx = d.ns[0];
	int ny = d.ns[1];
	int nz = d.ns[2];

	double h_x = d.spacings[0];
	double h_y = d.spacings[1];
	double h_z = d.spacings[2];

	LocalData<3>       f_data = f->getLocalData(d.local_index);
	const LocalData<3> u_data = u->getLocalData(d.local_index);

	double center, north, east, south, west, bottom, top;

	// derive in x direction
	// west
	if (d.hasNbr(Side<3>::west)) {
		const LocalData<2> boundary_data = gamma->getLocalData(d.getIfaceLocalIndex(Side<3>::west));
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                = boundary_data[{yi, zi}];
				center              = u_data[{0, yi, zi}];
				east                = u_data[{1, yi, zi}];
				f_data[{0, yi, zi}] = (2 * west - 3 * center + east) / (h_x * h_x);
			}
		}
	} else if (d.isNeumann(Side<3>::west)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				center              = u_data[{0, yi, zi}];
				east                = u_data[{1, yi, zi}];
				f_data[{0, yi, zi}] = (-center + east) / (h_x * h_x);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				center              = u_data[{0, yi, zi}];
				east                = u_data[{1, yi, zi}];
				f_data[{0, yi, zi}] = (-3 * center + east) / (h_x * h_x);
			}
		}
	}
	// middle
	for (int zi = 0; zi < nz; zi++) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 1; xi < nx - 1; xi++) {
				west   = u_data[{xi - 1, yi, zi}];
				center = u_data[{xi, yi, zi}];
				east   = u_data[{xi + 1, yi, zi}];

				f_data[{xi, yi, zi}] = (west - 2 * center + east) / (h_x * h_x);
			}
		}
	}
	// east
	if (d.hasNbr(Side<3>::east)) {
		const LocalData<2> boundary_data = gamma->getLocalData(d.getIfaceLocalIndex(Side<3>::east));
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                     = u_data[{nx - 2, yi, zi}];
				center                   = u_data[{nx - 1, yi, zi}];
				east                     = boundary_data[{yi, zi}];
				f_data[{nx - 1, yi, zi}] = (west - 3 * center + 2 * east) / (h_x * h_x);
			}
		}
	} else if (d.isNeumann(Side<3>::east)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                     = u_data[{nx - 2, yi, zi}];
				center                   = u_data[{nx - 1, yi, zi}];
				f_data[{nx - 1, yi, zi}] = (west - center) / (h_x * h_x);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                     = u_data[{nx - 2, yi, zi}];
				center                   = u_data[{nx - 1, yi, zi}];
				f_data[{nx - 1, yi, zi}] = (west - 3 * center) / (h_x * h_x);
			}
		}
	}

	// derive in y direction
	// south
	if (d.hasNbr(Side<3>::south)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(d.getIfaceLocalIndex(Side<3>::south));
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = boundary_data[{xi, zi}];
				center = u_data[{xi, 0, zi}];
				north  = u_data[{xi, 1, zi}];
				f_data[{xi, 0, zi}] += (2 * south - 3 * center + north) / (h_y * h_y);
			}
		}
	} else if (d.isNeumann(Side<3>::south)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u_data[{xi, 0, zi}];
				north  = u_data[{xi, 1, zi}];
				f_data[{xi, 0, zi}] += (-center + north) / (h_y * h_y);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u_data[{xi, 0, zi}];
				north  = u_data[{xi, 1, zi}];
				f_data[{xi, 0, zi}] += (-3 * center + north) / (h_y * h_y);
			}
		}
	}
	// middle
	for (int zi = 0; zi < nz; zi++) {
		for (int yi = 1; yi < ny - 1; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u_data[{xi, yi - 1, zi}];
				center = u_data[{xi, yi, zi}];
				north  = u_data[{xi, yi + 1, zi}];

				f_data[{xi, yi, zi}] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
	}
	// north
	if (d.hasNbr(Side<3>::north)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(d.getIfaceLocalIndex(Side<3>::north));
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u_data[{xi, ny - 2, zi}];
				center = u_data[{xi, ny - 1, zi}];
				north  = boundary_data[{xi, zi}];
				f_data[{xi, ny - 1, zi}] += (south - 3 * center + 2 * north) / (h_y * h_y);
			}
		}
	} else if (d.isNeumann(Side<3>::north)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u_data[{xi, ny - 2, zi}];
				center = u_data[{xi, ny - 1, zi}];
				f_data[{xi, ny - 1, zi}] += (south - center) / (h_y * h_y);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u_data[{xi, ny - 2, zi}];
				center = u_data[{xi, ny - 1, zi}];
				f_data[{xi, ny - 1, zi}] += (south - 3 * center) / (h_y * h_y);
			}
		}
	}

	// derive in z direction
	// bottom
	if (d.hasNbr(Side<3>::bottom)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(d.getIfaceLocalIndex(Side<3>::bottom));
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = boundary_data[{xi, yi}];
				center = u_data[{xi, yi, 0}];
				top    = u_data[{xi, yi, 1}];
				f_data[{xi, yi, 0}] += (2 * bottom - 3 * center + top) / (h_z * h_z);
			}
		}
	} else if (d.isNeumann(Side<3>::bottom)) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u_data[{xi, yi, 0}];
				top    = u_data[{xi, yi, 1}];
				f_data[{xi, yi, 0}] += (-center + top) / (h_z * h_z);
			}
		}
	} else {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u_data[{xi, yi, 0}];
				top    = u_data[{xi, yi, 1}];
				f_data[{xi, yi, 0}] += (-3 * center + top) / (h_z * h_z);
			}
		}
	}
	// middle
	for (int zi = 1; zi < nz - 1; zi++) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u_data[{xi, yi, zi - 1}];
				center = u_data[{xi, yi, zi}];
				top    = u_data[{xi, yi, zi + 1}];

				f_data[{xi, yi, zi}] += (bottom - 2 * center + top) / (h_z * h_z);
			}
		}
	}
	// top
	if (d.hasNbr(Side<3>::top)) {
		const LocalData<2> boundary_data = gamma->getLocalData(d.getIfaceLocalIndex(Side<3>::top));
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u_data[{xi, yi, nz - 2}];
				center = u_data[{xi, yi, nz - 1}];
				top    = boundary_data[{xi, yi}];
				f_data[{xi, yi, nz - 1}] += (bottom - 3 * center + 2 * top) / (h_z * h_z);
			}
		}
	} else if (d.isNeumann(Side<3>::top)) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u_data[{xi, yi, nz - 2}];
				center = u_data[{xi, yi, nz - 1}];
				f_data[{xi, yi, nz - 1}] += (bottom - center) / (h_z * h_z);
			}
		}
	} else {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u_data[{xi, yi, nz - 2}];
				center = u_data[{xi, yi, nz - 1}];
				f_data[{xi, yi, nz - 1}] += (bottom - 3 * center) / (h_z * h_z);
			}
		}
	}
}
