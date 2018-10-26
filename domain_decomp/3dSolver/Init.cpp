#include "Init.h"
#include <algorithm>
using namespace std;
void getXYZ(const Domain<3> &d, const int &xi, const int &yi, const int &zi, double &x, double &y,
            double &z)
{
	const int &n   = d.n;
	double     h_x = d.lengths[0] / n;
	double     h_y = d.lengths[1] / n;
	double     h_z = d.lengths[2] / n;
	if (xi == -1) {
		x = d.starts[0];
	} else if (xi == n) {
		x = d.starts[0] + d.lengths[0];
	} else {
		x = d.starts[0] + h_x / 2.0 + d.lengths[0] * xi / n;
	}
	if (yi == -1) {
		y = d.starts[1];
	} else if (yi == n) {
		y = d.starts[1] + d.lengths[1];
	} else {
		y = d.starts[1] + h_y / 2.0 + d.lengths[1] * yi / n;
	}
	if (zi == -1) {
		z = d.starts[2];
	} else if (zi == n) {
		z = d.starts[2] + d.lengths[2];
	} else {
		z = d.starts[2] + h_z / 2.0 + d.lengths[2] * zi / n;
	}
}
inline int index(Domain<3> &d, const int &xi, const int &yi, const int &zi)
{
	const int &n = d.n;
	return xi + yi * n + zi * n * n;
}
void Init::initNeumann(DomainCollection &dc, int n, Vec f, Vec exact,
                       function<double(double, double, double)> ffun,
                       function<double(double, double, double)> efun,
                       function<double(double, double, double)> nfunx,
                       function<double(double, double, double)> nfuny,
                       function<double(double, double, double)> nfunz)
{
	n = dc.getN();
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &p : dc.domains) {
		Domain<3> &d = *p.second;

		double *f_vals     = f_ptr + d.id_local * n * n * n;
		double *exact_vals = exact_ptr + d.id_local * n * n * n;

		// Generate RHS vector
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, zi, x, y, z);
					f_vals[index(d, xi, yi, zi)]     = ffun(x, y, z);
					exact_vals[index(d, xi, yi, zi)] = efun(x, y, z);
				}
			}
		}
		// apply boundaries
		double h_x = d.lengths[0] / n;
		double h_y = d.lengths[1] / n;
		double h_z = d.lengths[2] / n;
		// west
		if (!d.hasNbr(Side<3>::west)) {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					double x, y, z;
					getXYZ(d, -1, yi, zi, x, y, z);
					f_vals[index(d, 0, yi, zi)] += nfunx(x, y, z) / h_x;
				}
			}
		}
		// east
		if (!d.hasNbr(Side<3>::east)) {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					double x, y, z;
					getXYZ(d, n, yi, zi, x, y, z);
					f_vals[index(d, n - 1, yi, zi)] -= nfunx(x, y, z) / h_x;
				}
			}
		}
		// south
		if (!d.hasNbr(Side<3>::south)) {
			for (int zi = 0; zi < n; zi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, -1, zi, x, y, z);
					f_vals[index(d, xi, 0, zi)] += nfuny(x, y, z) / h_y;
				}
			}
		}
		// north
		if (!d.hasNbr(Side<3>::north)) {
			for (int zi = 0; zi < n; zi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, n, zi, x, y, z);
					f_vals[index(d, xi, n - 1, zi)] -= nfuny(x, y, z) / h_y;
				}
			}
		}
		// bottom
		if (!d.hasNbr(Side<3>::bottom)) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, -1, x, y, z);
					f_vals[index(d, xi, yi, 0)] += nfunz(x, y, z) / h_z;
				}
			}
		}
		// top
		if (!d.hasNbr(Side<3>::top)) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, n, x, y, z);
					f_vals[index(d, xi, yi, n - 1)] -= nfunz(x, y, z) / h_z;
				}
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::initDirichlet(DomainCollection &dc, int n, Vec f, Vec exact,
                         function<double(double, double, double)> ffun,
                         function<double(double, double, double)> efun)
{
	n = dc.getN();
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &p : dc.domains) {
		Domain<3> &d = *p.second;

		double *f_vals     = f_ptr + d.id_local * n * n * n;
		double *exact_vals = exact_ptr + d.id_local * n * n * n;
		// Generate RHS vector
		double h_x = d.lengths[0] / n;
		double h_y = d.lengths[1] / n;
		double h_z = d.lengths[2] / n;
		h_x *= h_x;
		h_y *= h_y;
		h_z *= h_z;
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, zi, x, y, z);
					f_vals[index(d, xi, yi, zi)]     = ffun(x, y, z);
					exact_vals[index(d, xi, yi, zi)] = efun(x, y, z);
				}
			}
		}
		// apply boundaries
		// west
		if (!d.hasNbr(Side<3>::west)) {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					double x, y, z;
					getXYZ(d, -1, yi, zi, x, y, z);
					f_vals[index(d, 0, yi, zi)] -= 2 * efun(x, y, z) / h_x;
				}
			}
		}
		// east
		if (!d.hasNbr(Side<3>::east)) {
			for (int zi = 0; zi < n; zi++) {
				for (int yi = 0; yi < n; yi++) {
					double x, y, z;
					getXYZ(d, n, yi, zi, x, y, z);
					f_vals[index(d, n - 1, yi, zi)] -= 2 * efun(x, y, z) / h_x;
				}
			}
		}
		// south
		if (!d.hasNbr(Side<3>::south)) {
			for (int zi = 0; zi < n; zi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, -1, zi, x, y, z);
					f_vals[index(d, xi, 0, zi)] -= 2 * efun(x, y, z) / h_y;
				}
			}
		}
		// north
		if (!d.hasNbr(Side<3>::north)) {
			for (int zi = 0; zi < n; zi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, n, zi, x, y, z);
					f_vals[index(d, xi, n - 1, zi)] -= 2 * efun(x, y, z) / h_y;
				}
			}
		}
		// bottom
		if (!d.hasNbr(Side<3>::bottom)) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, -1, x, y, z);
					f_vals[index(d, xi, yi, 0)] -= 2 * efun(x, y, z) / h_z;
				}
			}
		}
		// top
		if (!d.hasNbr(Side<3>::top)) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					double x, y, z;
					getXYZ(d, xi, yi, n, x, y, z);
					f_vals[index(d, xi, yi, n - 1)] -= 2 * efun(x, y, z) / h_z;
				}
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
