#include "Init.h"
using namespace std;
void Init::initNeumann(DomainSignatureCollection &dsc, int n, double *f_vec, double *exact_vec,
                       function<double(double, double)> ffun, function<double(double, double)> efun,
                       function<double(double, double)> nfunx,
                       function<double(double, double)> nfuny)
{
	for (auto &p : dsc.domains) {
		DomainSignature &d = p.second;

		// Generate RHS vector
		double *f     = f_vec + d.id * n * n;
		double *exact = exact_vec + d.id * n * n;

		double h_x = d.x_length / n;
		double h_y = d.y_length / n;
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x           = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = efun(x, y);
			}
		}
		// apply boundaries
		// west
		if (!d.hasNbr(Side::west)) {
			for (int yi = 0; yi < n; yi++) {
				double y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n] += nfunx(d.x_start, y) / h_x;
			}
		}
		// east
		if (!d.hasNbr(Side::east)) {
			for (int yi = 0; yi < n; yi++) {
				double y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + n - 1] -= nfunx(d.x_start + d.x_length, y) / h_x;
			}
		}
		// south
		if (!d.hasNbr(Side::south)) {
			for (int xi = 0; xi < n; xi++) {
				double x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				f[xi] += nfuny(x, d.y_start) / h_y;
			}
		}
		// north
		if (!d.hasNbr(Side::north)) {
			for (int xi = 0; xi < n; xi++) {
				double x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				f[n * (n - 1) + xi] -= nfuny(x, d.y_start+d.y_length) / h_y;
			}
		}
	}
}
void Init::initDirichlet(DomainSignatureCollection &dsc, int n, double *f_vec, double *exact_vec,
                       function<double(double, double)> ffun, function<double(double, double)> efun)
{
    n = dsc.n;
	for (auto &p : dsc.domains) {
		DomainSignature &d = p.second;

		// Generate RHS vector
		double *f     = f_vec + d.id * n * n;
		double *exact = exact_vec + d.id * n * n;

		double h_x = d.x_length / n;
		double h_y = d.y_length / n;
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x           = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = efun(x, y);
			}
		}
		// apply boundaries
		// west
		if (!d.hasNbr(Side::west)) {
			for (int yi = 0; yi < n; yi++) {
				double y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n] -= efun(d.x_start, y)*2 / (h_x*h_x);
			}
		}
		// east
		if (!d.hasNbr(Side::east)) {
			for (int yi = 0; yi < n; yi++) {
				double y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + n - 1] -= efun(d.x_start + d.x_length, y) *2 / (h_x*h_x);
			}
		}
		// south
		if (!d.hasNbr(Side::south)) {
			for (int xi = 0; xi < n; xi++) {
				double x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				f[xi] -= efun(x, d.y_start)*2 / (h_y*h_y);
			}
		}
		// north
		if (!d.hasNbr(Side::north)) {
			for (int xi = 0; xi < n; xi++) {
				double x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				f[n * (n - 1) + xi] -= efun(x, d.y_start+d.y_length)*2 / (h_y*h_y);
			}
		}
	}
}
