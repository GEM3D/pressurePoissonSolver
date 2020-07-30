#include "SevenPtPatchOperator.h"
void SevenPtPatchOperator::applyWithInterface(SchurInfo<3> &sinfo, const LocalData<3> u,
                                              std::shared_ptr<const Vector<2>> gamma,
                                              LocalData<3>                     f)

{
	int nx = sinfo.pinfo->ns[0];
	int ny = sinfo.pinfo->ns[1];
	int nz = sinfo.pinfo->ns[2];

	double h_x = sinfo.pinfo->spacings[0];
	double h_y = sinfo.pinfo->spacings[1];
	double h_z = sinfo.pinfo->spacings[2];

	double center, north, east, south, west, bottom, top;

	// derive in x direction
	// west
	if (sinfo.pinfo->hasNbr(Side<3>::west)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<3>::west));
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west           = boundary_data[{yi, zi}];
				center         = u[{0, yi, zi}];
				east           = u[{1, yi, zi}];
				f[{0, yi, zi}] = (2 * west - 3 * center + east) / (h_x * h_x);
			}
		}
	} else if (sinfo.pinfo->isNeumann(Side<3>::west)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				center         = u[{0, yi, zi}];
				east           = u[{1, yi, zi}];
				f[{0, yi, zi}] = (-center + east) / (h_x * h_x);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				center         = u[{0, yi, zi}];
				east           = u[{1, yi, zi}];
				f[{0, yi, zi}] = (-3 * center + east) / (h_x * h_x);
			}
		}
	}
	// middle
	for (int zi = 0; zi < nz; zi++) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 1; xi < nx - 1; xi++) {
				west   = u[{xi - 1, yi, zi}];
				center = u[{xi, yi, zi}];
				east   = u[{xi + 1, yi, zi}];

				f[{xi, yi, zi}] = (west - 2 * center + east) / (h_x * h_x);
			}
		}
	}
	// east
	if (sinfo.pinfo->hasNbr(Side<3>::east)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<3>::east));
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                = u[{nx - 2, yi, zi}];
				center              = u[{nx - 1, yi, zi}];
				east                = boundary_data[{yi, zi}];
				f[{nx - 1, yi, zi}] = (west - 3 * center + 2 * east) / (h_x * h_x);
			}
		}
	} else if (sinfo.pinfo->isNeumann(Side<3>::east)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                = u[{nx - 2, yi, zi}];
				center              = u[{nx - 1, yi, zi}];
				f[{nx - 1, yi, zi}] = (west - center) / (h_x * h_x);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                = u[{nx - 2, yi, zi}];
				center              = u[{nx - 1, yi, zi}];
				f[{nx - 1, yi, zi}] = (west - 3 * center) / (h_x * h_x);
			}
		}
	}

	// derive in y direction
	// south
	if (sinfo.pinfo->hasNbr(Side<3>::south)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<3>::south));
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = boundary_data[{xi, zi}];
				center = u[{xi, 0, zi}];
				north  = u[{xi, 1, zi}];
				f[{xi, 0, zi}] += (2 * south - 3 * center + north) / (h_y * h_y);
			}
		}
	} else if (sinfo.pinfo->isNeumann(Side<3>::south)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, 0, zi}];
				north  = u[{xi, 1, zi}];
				f[{xi, 0, zi}] += (-center + north) / (h_y * h_y);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, 0, zi}];
				north  = u[{xi, 1, zi}];
				f[{xi, 0, zi}] += (-3 * center + north) / (h_y * h_y);
			}
		}
	}
	// middle
	for (int zi = 0; zi < nz; zi++) {
		for (int yi = 1; yi < ny - 1; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, yi - 1, zi}];
				center = u[{xi, yi, zi}];
				north  = u[{xi, yi + 1, zi}];

				f[{xi, yi, zi}] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
	}
	// north
	if (sinfo.pinfo->hasNbr(Side<3>::north)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<3>::north));
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, ny - 2, zi}];
				center = u[{xi, ny - 1, zi}];
				north  = boundary_data[{xi, zi}];
				f[{xi, ny - 1, zi}] += (south - 3 * center + 2 * north) / (h_y * h_y);
			}
		}
	} else if (sinfo.pinfo->isNeumann(Side<3>::north)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, ny - 2, zi}];
				center = u[{xi, ny - 1, zi}];
				f[{xi, ny - 1, zi}] += (south - center) / (h_y * h_y);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, ny - 2, zi}];
				center = u[{xi, ny - 1, zi}];
				f[{xi, ny - 1, zi}] += (south - 3 * center) / (h_y * h_y);
			}
		}
	}

	// derive in z direction
	// bottom
	if (sinfo.pinfo->hasNbr(Side<3>::bottom)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<3>::bottom));
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = boundary_data[{xi, yi}];
				center = u[{xi, yi, 0}];
				top    = u[{xi, yi, 1}];
				f[{xi, yi, 0}] += (2 * bottom - 3 * center + top) / (h_z * h_z);
			}
		}
	} else if (sinfo.pinfo->isNeumann(Side<3>::bottom)) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, yi, 0}];
				top    = u[{xi, yi, 1}];
				f[{xi, yi, 0}] += (-center + top) / (h_z * h_z);
			}
		}
	} else {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, yi, 0}];
				top    = u[{xi, yi, 1}];
				f[{xi, yi, 0}] += (-3 * center + top) / (h_z * h_z);
			}
		}
	}
	// middle
	for (int zi = 1; zi < nz - 1; zi++) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u[{xi, yi, zi - 1}];
				center = u[{xi, yi, zi}];
				top    = u[{xi, yi, zi + 1}];

				f[{xi, yi, zi}] += (bottom - 2 * center + top) / (h_z * h_z);
			}
		}
	}
	// top
	if (sinfo.pinfo->hasNbr(Side<3>::top)) {
		const LocalData<2> boundary_data
		= gamma->getLocalData(sinfo.getIfaceLocalIndex(Side<3>::top));
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u[{xi, yi, nz - 2}];
				center = u[{xi, yi, nz - 1}];
				top    = boundary_data[{xi, yi}];
				f[{xi, yi, nz - 1}] += (bottom - 3 * center + 2 * top) / (h_z * h_z);
			}
		}
	} else if (sinfo.pinfo->isNeumann(Side<3>::top)) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u[{xi, yi, nz - 2}];
				center = u[{xi, yi, nz - 1}];
				f[{xi, yi, nz - 1}] += (bottom - center) / (h_z * h_z);
			}
		}
	} else {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u[{xi, yi, nz - 2}];
				center = u[{xi, yi, nz - 1}];
				f[{xi, yi, nz - 1}] += (bottom - 3 * center) / (h_z * h_z);
			}
		}
	}
}
void SevenPtPatchOperator::addInterfaceToRHS(SchurInfo<3> &                   sinfo,
                                             std::shared_ptr<const Vector<2>> gamma, LocalData<3> f)
{
	for (Side<3> s : Side<3>::getValues()) {
		if (sinfo.pinfo->hasNbr(s)) {
			const LocalData<2> gamma_view = gamma->getLocalData(sinfo.getIfaceLocalIndex(s));
			LocalData<2>       slice      = f.getSliceOnSide(s);
			double             h2         = pow(sinfo.pinfo->spacings[s.axis()], 2);
			nested_loop<2>(
			gamma_view.getStart(), gamma_view.getEnd(),
			[&](std::array<int, 2> coord) { slice[coord] -= 2.0 / h2 * gamma_view[coord]; });
		}
	}
}
void SevenPtPatchOperator::apply(const SchurInfo<3> &sinfo, const LocalData<3> u, LocalData<3> f)

{
	int nx = sinfo.pinfo->ns[0];
	int ny = sinfo.pinfo->ns[1];
	int nz = sinfo.pinfo->ns[2];

	double h_x = sinfo.pinfo->spacings[0];
	double h_y = sinfo.pinfo->spacings[1];
	double h_z = sinfo.pinfo->spacings[2];

	double center, north, east, south, west, bottom, top;

	// derive in x direction
	// west
	if (sinfo.pinfo->isNeumann(Side<3>::west)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				center         = u[{0, yi, zi}];
				east           = u[{1, yi, zi}];
				f[{0, yi, zi}] = (-center + east) / (h_x * h_x);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				center         = u[{0, yi, zi}];
				east           = u[{1, yi, zi}];
				f[{0, yi, zi}] = (-3 * center + east) / (h_x * h_x);
			}
		}
	}
	// middle
	for (int zi = 0; zi < nz; zi++) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 1; xi < nx - 1; xi++) {
				west   = u[{xi - 1, yi, zi}];
				center = u[{xi, yi, zi}];
				east   = u[{xi + 1, yi, zi}];

				f[{xi, yi, zi}] = (west - 2 * center + east) / (h_x * h_x);
			}
		}
	}
	// east
	if (sinfo.pinfo->isNeumann(Side<3>::east)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                = u[{nx - 2, yi, zi}];
				center              = u[{nx - 1, yi, zi}];
				f[{nx - 1, yi, zi}] = (west - center) / (h_x * h_x);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				west                = u[{nx - 2, yi, zi}];
				center              = u[{nx - 1, yi, zi}];
				f[{nx - 1, yi, zi}] = (west - 3 * center) / (h_x * h_x);
			}
		}
	}

	// derive in y direction
	// south
	if (sinfo.pinfo->isNeumann(Side<3>::south)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, 0, zi}];
				north  = u[{xi, 1, zi}];
				f[{xi, 0, zi}] += (-center + north) / (h_y * h_y);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, 0, zi}];
				north  = u[{xi, 1, zi}];
				f[{xi, 0, zi}] += (-3 * center + north) / (h_y * h_y);
			}
		}
	}
	// middle
	for (int zi = 0; zi < nz; zi++) {
		for (int yi = 1; yi < ny - 1; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, yi - 1, zi}];
				center = u[{xi, yi, zi}];
				north  = u[{xi, yi + 1, zi}];

				f[{xi, yi, zi}] += (south - 2 * center + north) / (h_y * h_y);
			}
		}
	}
	// north
	if (sinfo.pinfo->isNeumann(Side<3>::north)) {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, ny - 2, zi}];
				center = u[{xi, ny - 1, zi}];
				f[{xi, ny - 1, zi}] += (south - center) / (h_y * h_y);
			}
		}
	} else {
		for (int zi = 0; zi < nz; zi++) {
			for (int xi = 0; xi < nx; xi++) {
				south  = u[{xi, ny - 2, zi}];
				center = u[{xi, ny - 1, zi}];
				f[{xi, ny - 1, zi}] += (south - 3 * center) / (h_y * h_y);
			}
		}
	}

	// derive in z direction
	// bottom
	if (sinfo.pinfo->isNeumann(Side<3>::bottom)) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, yi, 0}];
				top    = u[{xi, yi, 1}];
				f[{xi, yi, 0}] += (-center + top) / (h_z * h_z);
			}
		}
	} else {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				center = u[{xi, yi, 0}];
				top    = u[{xi, yi, 1}];
				f[{xi, yi, 0}] += (-3 * center + top) / (h_z * h_z);
			}
		}
	}
	// middle
	for (int zi = 1; zi < nz - 1; zi++) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u[{xi, yi, zi - 1}];
				center = u[{xi, yi, zi}];
				top    = u[{xi, yi, zi + 1}];

				f[{xi, yi, zi}] += (bottom - 2 * center + top) / (h_z * h_z);
			}
		}
	}
	// top
	if (sinfo.pinfo->isNeumann(Side<3>::top)) {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u[{xi, yi, nz - 2}];
				center = u[{xi, yi, nz - 1}];
				f[{xi, yi, nz - 1}] += (bottom - center) / (h_z * h_z);
			}
		}
	} else {
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				bottom = u[{xi, yi, nz - 2}];
				center = u[{xi, yi, nz - 1}];
				f[{xi, yi, nz - 1}] += (bottom - 3 * center) / (h_z * h_z);
			}
		}
	}
}
