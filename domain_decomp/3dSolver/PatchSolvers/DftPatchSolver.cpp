#include "DftPatchSolver.h"
#include "Utils.h"
using namespace std;
using namespace Utils;

extern "C" void dgemv_(char &, int &, int &, double &, double *, int &, double *, int &, double &,
                       double *, int &);
inline int      index(const int &n, const int &xi, const int &yi, const int &zi)
{
	return xi + yi * n + zi * n * n;
}
DftPatchSolver::DftPatchSolver(DomainCollection &dc, double lambda)
{
	n            = dc.getN();
	this->lambda = lambda;
}
void DftPatchSolver::addDomain(SchurDomain<3> &d)
{
	if (!initialized) {
		initialized = true;
		f_copy.resize(n * n * n);
		tmp.resize(n * n * n);
	}

	if (!plan1.count(d)) {
		DftType x_transform     = DftType::DST_II;
		DftType x_transform_inv = DftType::DST_III;
		DftType y_transform     = DftType::DST_II;
		DftType y_transform_inv = DftType::DST_III;
		DftType z_transform     = DftType::DST_II;
		DftType z_transform_inv = DftType::DST_III;
		if (d.neumann.to_ulong()) {
			// x direction
			if (d.isNeumann(Side<3>::east) && d.isNeumann(Side<3>::west)) {
				x_transform     = DftType::DCT_II;
				x_transform_inv = DftType::DCT_III;
			} else if (d.isNeumann(Side<3>::west)) {
				x_transform     = DftType::DCT_IV;
				x_transform_inv = DftType::DCT_IV;
			} else if (d.isNeumann(Side<3>::east)) {
				x_transform     = DftType::DST_IV;
				x_transform_inv = DftType::DST_IV;
			}
			// y direction
			if (d.isNeumann(Side<3>::north) && d.isNeumann(Side<3>::south)) {
				y_transform     = DftType::DCT_II;
				y_transform_inv = DftType::DCT_III;
			} else if (d.isNeumann(Side<3>::south)) {
				y_transform     = DftType::DCT_IV;
				y_transform_inv = DftType::DCT_IV;
			} else if (d.isNeumann(Side<3>::north)) {
				y_transform     = DftType::DST_IV;
				y_transform_inv = DftType::DST_IV;
			}
			// z direction
			if (d.isNeumann(Side<3>::bottom) && d.isNeumann(Side<3>::top)) {
				z_transform     = DftType::DCT_II;
				z_transform_inv = DftType::DCT_III;
			} else if (d.isNeumann(Side<3>::bottom)) {
				z_transform     = DftType::DCT_IV;
				z_transform_inv = DftType::DCT_IV;
			} else if (d.isNeumann(Side<3>::top)) {
				z_transform     = DftType::DST_IV;
				z_transform_inv = DftType::DST_IV;
			}
		}

		plan1[d] = plan(x_transform, y_transform, z_transform);
		plan2[d] = plan(x_transform_inv, y_transform_inv, z_transform_inv);
	}

	double h_x = d.domain.lengths[0] / n;
	double h_y = d.domain.lengths[1] / n;
	double h_z = d.domain.lengths[2] / n;
	if (!eigen_vals.count(d)) {
		valarray<double> &eigen_val = eigen_vals[d];
		eigen_val.resize(n * n * n);
		// create eigen_val vector
		// z direction
		if (d.isNeumann(Side<3>::bottom) && d.isNeumann(Side<3>::top)) {
			for (int zi = 0; zi < n; zi++) {
				eigen_val[slice(zi * n * n, n * n, 1)]
				= -4 / (h_z * h_z) * pow(sin(zi * M_PI / (2 * n)), 2);
			}
		} else if (d.isNeumann(Side<3>::bottom) || d.isNeumann(Side<3>::top)) {
			for (int zi = 0; zi < n; zi++) {
				eigen_val[slice(zi * n * n, n * n, 1)]
				= -4 / (h_z * h_z) * pow(sin((zi + 0.5) * M_PI / (2 * n)), 2);
			}
		} else {
			for (int zi = 0; zi < n; zi++) {
				eigen_val[slice(zi * n * n, n * n, 1)]
				= -4 / (h_z * h_z) * pow(sin((zi + 1) * M_PI / (2 * n)), 2);
			}
		}

		valarray<double> ones(n * n);
		ones = 1;

		// y direction
		valarray<size_t> sizes = {(size_t) n, (size_t) n};
		valarray<size_t> strides(2);
		strides[0] = n * n;
		strides[1] = 1;
		if (d.isNeumann(Side<3>::south) && d.isNeumann(Side<3>::north)) {
			for (int yi = 0; yi < n; yi++) {
				eigen_val[gslice(yi * n, sizes, strides)]
				-= 4 / (h_y * h_y) * pow(sin(yi * M_PI / (2 * n)), 2) * ones;
			}
		} else if (d.isNeumann(Side<3>::south) || d.isNeumann(Side<3>::north)) {
			for (int yi = 0; yi < n; yi++) {
				eigen_val[gslice(yi * n, sizes, strides)]
				-= 4 / (h_y * h_y) * pow(sin((yi + 0.5) * M_PI / (2 * n)), 2) * ones;
			}
		} else {
			for (int yi = 0; yi < n; yi++) {
				eigen_val[gslice(yi * n, sizes, strides)]
				-= 4 / (h_y * h_y) * pow(sin((yi + 1) * M_PI / (2 * n)), 2) * ones;
			}
		}

		// x direction
		strides[0] = n * n;
		strides[1] = n;
		if (d.isNeumann(Side<3>::west) && d.isNeumann(Side<3>::east)) {
			for (int xi = 0; xi < n; xi++) {
				eigen_val[gslice(xi, sizes, strides)]
				-= 4 / (h_x * h_x) * pow(sin(xi * M_PI / (2 * n)), 2) * ones;
			}
		} else if (d.isNeumann(Side<3>::west) || d.isNeumann(Side<3>::east)) {
			for (int xi = 0; xi < n; xi++) {
				eigen_val[gslice(xi, sizes, strides)]
				-= 4 / (h_x * h_x) * pow(sin((xi + 0.5) * M_PI / (2 * n)), 2) * ones;
			}
		} else {
			for (int xi = 0; xi < n; xi++) {
				eigen_val[gslice(xi, sizes, strides)]
				-= 4 / (h_x * h_x) * pow(sin((xi + 1) * M_PI / (2 * n)), 2) * ones;
			}
		}
		eigen_val += lambda;
	}
}
void DftPatchSolver::solve(SchurDomain<3> &d, const Vec f, Vec u, const Vec gamma)
{
	double h_x        = d.domain.lengths[0] / n;
	double h_y        = d.domain.lengths[1] / n;
	double h_z        = d.domain.lengths[2] / n;
	auto   getSpacing = [=](Side<3> s) {
        double retval = 0;
        switch (s.toInt()) {
            case Side<3>::east:
            case Side<3>::west:
                retval = h_x;
                break;
            case Side<3>::south:
            case Side<3>::north:
                retval = h_y;
                break;
            case Side<3>::bottom:
            case Side<3>::top:
                retval = h_z;
        }
        return retval;
	};

	const double *f_view, *gamma_view;
	VecGetArrayRead(f, &f_view);
	VecGetArrayRead(gamma, &gamma_view);

	int start = d.local_index * n * n * n;
	for (int i = 0; i < n * n * n; i++) {
		f_copy[i] = f_view[start + i];
	}

	for (Side<3> s : Side<3>::getValues()) {
		if (d.hasNbr(s)) {
			int    idx = n * n * d.getIfaceLocalIndex(s);
			Slice  sl  = getSlice(&f_copy[0], n, s);
			double h2  = pow(getSpacing(s), 2);
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					sl(xi, yi) -= 2.0 / h2 * gamma_view[idx + xi + yi * n];
				}
			}
		}
	}

	double *u_view;
	VecGetArray(u, &u_view);
	double *u_local_view = u_view + start;
	tmp                  = 0;

	execute_plan(plan1[d], &f_copy[0], &tmp[0], false);

	tmp /= eigen_vals[d];

	if (d.neumann.all()) { tmp[0] = 0; }

	execute_plan(plan2[d], &tmp[0], u_local_view, true);

	double scale = 8.0 / pow(n, 3);
	int    num   = pow(n, 3);
	for (int i = 0; i < num; i++) {
		u_local_view[i] *= scale;
	}
	VecRestoreArray(u, &u_view);
	VecRestoreArrayRead(f, &f_view);
	VecRestoreArrayRead(gamma, &gamma_view);
}
std::array<std::shared_ptr<valarray<double>>, 3>
DftPatchSolver::plan(DftType x_type, DftType y_type, DftType z_type)
{
	array<shared_ptr<valarray<double>>, 3> retval;
	retval[0] = getTransformArray(x_type);
	retval[1] = getTransformArray(y_type);
	retval[2] = getTransformArray(z_type);
	return retval;
}
std::shared_ptr<valarray<double>> DftPatchSolver::getTransformArray(DftType type)
{
	int idx = static_cast<int>(type);

	shared_ptr<valarray<double>> matrix_ptr;

	if (transforms[idx] == nullptr) {
		matrix_ptr.reset(new valarray<double>(n * n));
		valarray<double> &matrix = *matrix_ptr;
		switch (type) {
			case DftType::DCT_II:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * (i * (j + 0.5)));
					}
				}
				break;
			case DftType::DCT_III:
				for (int i = 0; i < n; i++) {
					matrix[i * n] = 0.5;
				}
				for (int j = 1; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * ((i + 0.5) * j));
					}
				}
				break;
			case DftType::DCT_IV:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = cos(M_PI / n * ((i + 0.5) * (j + 0.5)));
					}
				}
				break;
			case DftType::DST_II:
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 1) * (j + 0.5)));
					}
				}
				break;
			case DftType::DST_III:
				for (int i = 0; i < n; i += 2) {
					matrix[i * n + n - 1] = 0.5;
				}
				for (int i = 1; i < n; i += 2) {
					matrix[i * n + n - 1] = -0.5;
				}
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n - 1; j++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 0.5) * (j + 1)));
					}
				}
				break;
			case DftType::DST_IV:
				for (int j = 0; j < n; j++) {
					for (int i = 0; i < n; i++) {
						matrix[i * n + j] = sin(M_PI / n * ((i + 0.5) * (j + 0.5)));
					}
				}
				break;
		}
	} else {
		matrix_ptr = transforms[idx];
	}
	return matrix_ptr;
}
void DftPatchSolver::execute_plan(std::array<std::shared_ptr<std::valarray<double>>, 3> plan,
                                  double *in, double *out, const bool inverse)
{
	double *  prev_result = in;
	const int strides[3]  = {1, n, n * n};
	std::fill(out, out + n * n * n, 0);
	for (int dim = 0; dim < 3; dim++) {
		int other_strides[2];
		for (int i = 0; i < dim; i++) {
			other_strides[i] = strides[i];
		}
		int dft_stride = strides[dim];
		for (int i = dim + 1; i < 3; i++) {
			other_strides[i - 1] = strides[i];
		}
		int x_stride = other_strides[0];
		int y_stride = other_strides[1];

		std::valarray<double> &matrix = *plan[dim];
		double *               new_result;
		if (dim != 2) {
			new_result = new double[n * n * n]();
		} else {
			new_result = out;
		}
		for (int y = 0; y < n; y++) {
			for (int x = 0; x < n; x++) {
#if 1
				char   T   = 'T';
				double one = 1;
				dgemv_(T, n, n, one, &matrix[0], n, &prev_result[y * y_stride + x * x_stride],
				       dft_stride, one, &new_result[y * y_stride + x * x_stride], dft_stride);
#else
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						new_result[y * y_stride + x * x_stride + i * dft_stride]
						+= matrix[i * n + j]
						   * prev_result[y * y_stride + x * x_stride + j * dft_stride];
					}
				}
#endif
			}
		}
		if (dim != 0) { delete[] prev_result; }
		prev_result = new_result;
	}
}
