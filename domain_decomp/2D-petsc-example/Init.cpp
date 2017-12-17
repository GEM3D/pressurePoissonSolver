#include "Init.h"
using namespace std;
void Init::initNeumann(DomainCollection &dc, int n, Vec f, Vec exact,
                       function<double(double, double)> ffun, function<double(double, double)> efun,
                       function<double(double, double)> nfunx,
                       function<double(double, double)> nfuny)
{
	n = dc.n;
	vector<int> inds(dc.domains.size() * n * n);
	iota(inds.begin(), inds.end(), 0);
	vector<double> f_vals(dc.domains.size() * n * n);
	vector<double> exact_vals(dc.domains.size() * n * n);
	for (auto &p : dc.domains) {
		Domain &d = p.second;

		// Generate RHS vector
		double h_x = d.x_length / n;
		double h_y = d.y_length / n;
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x                = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				double y                = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f_vals[yi * n + xi]     = ffun(x, y);
				exact_vals[yi * n + xi] = efun(x, y);
			}
		}
		// apply boundaries
		// west
		if (!d.hasNbr(Side::west)) {
			for (int yi = 0; yi < n; yi++) {
				double y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f_vals[yi * n] += nfunx(d.x_start, y) / h_x;
			}
		}
		// east
		if (!d.hasNbr(Side::east)) {
			for (int yi = 0; yi < n; yi++) {
				double y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f_vals[yi * n + n - 1] -= nfunx(d.x_start + d.x_length, y) / h_x;
			}
		}
		// south
		if (!d.hasNbr(Side::south)) {
			for (int xi = 0; xi < n; xi++) {
				double x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				f_vals[xi] += nfuny(x, d.y_start) / h_y;
			}
		}
		// north
		if (!d.hasNbr(Side::north)) {
			for (int xi = 0; xi < n; xi++) {
				double x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				f_vals[n * (n - 1) + xi] -= nfuny(x, d.y_start + d.y_length) / h_y;
			}
		}
	}
	VecSetValuesLocal(f, inds.size(), &inds[0], &f_vals[0], INSERT_VALUES);
	VecSetValuesLocal(exact, inds.size(), &inds[0], &exact_vals[0], INSERT_VALUES);
	VecAssemblyBegin(f);
	VecAssemblyBegin(exact);
	VecAssemblyEnd(f);
	VecAssemblyEnd(exact);
}
void Init::initDirichlet(DomainCollection &dc, int n, Vec f, Vec exact,
                         function<double(double, double)> ffun,
                         function<double(double, double)> efun)
{
	n = dc.n;
	double *f_ptr;
	VecGetArray(f, &f_ptr);
	double *exact_ptr;
	VecGetArray(exact, &exact_ptr);
	for (auto &p : dc.domains) {
		Domain &d = p.second;

		double *f_vals     = f_ptr + d.id_local * n * n;
		double *exact_vals = exact_ptr + d.id_local * n * n;
		// Generate RHS vector
		double h_x = d.x_length / n;
		double h_y = d.y_length / n;
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x                = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				double y                = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f_vals[yi * n + xi]     = ffun(x, y);
				exact_vals[yi * n + xi] = efun(x, y);
			}
		}
		// apply boundaries
		// west
		if (!d.hasNbr(Side::west)) {
			for (int yi = 0; yi < n; yi++) {
				double y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f_vals[yi * n] -= efun(d.x_start, y) * 2 / (h_x * h_x);
			}
		}
		// east
		if (!d.hasNbr(Side::east)) {
			for (int yi = 0; yi < n; yi++) {
				double y = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f_vals[yi * n + n - 1] -= efun(d.x_start + d.x_length, y) * 2 / (h_x * h_x);
			}
		}
		// south
		if (!d.hasNbr(Side::south)) {
			for (int xi = 0; xi < n; xi++) {
				double x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				f_vals[xi] -= efun(x, d.y_start) * 2 / (h_y * h_y);
			}
		}
		// north
		if (!d.hasNbr(Side::north)) {
			for (int xi = 0; xi < n; xi++) {
				double x = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				f_vals[n * (n - 1) + xi] -= efun(x, d.y_start + d.y_length) * 2 / (h_y * h_y);
			}
		}
	}
	VecRestoreArray(f, &f_ptr);
	VecRestoreArray(exact, &exact_ptr);
}
void Init::fillSolution(DomainCollection &dc, Vec u, function<double(double, double, double)> fun,
                        double time)
{
	double *vec;
	VecGetArray(u, &vec);
	int n = dc.n;
	for (auto &p : dc.domains) {
		Domain &d = p.second;

		double *f = vec + d.id_local * n * n;

		double h_x = d.x_length / n;
		double h_y = d.y_length / n;
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x       = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				double y       = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi] = fun(x, y, time);
			}
		}
	}
	VecRestoreArray(u, &vec);
}
