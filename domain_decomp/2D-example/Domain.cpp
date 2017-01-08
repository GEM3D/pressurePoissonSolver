#include "Domain.h"
#include <cmath>
#include <fftw3.h>
#include <iostream>
using namespace Eigen;
using namespace std;
Domain::Domain(MyArray grid, double h_x, double h_y)
{
	this->grid = grid;
	ny         = grid.rows();
	nx         = grid.cols();
	this->u    = MyArray(ny, nx);
	this->h_x  = h_x;
	this->h_y  = h_y;

	MyArray jx = RowVectorXd::LinSpaced(nx, 1, nx);
	jx = jx * (M_PI / (2 * nx));
	jx = jx.sin().square();
	jx = jx * (-4 / (h_x * h_x));

	MyArray jy = VectorXd::LinSpaced(ny, 1, ny);
	jy = jy * (M_PI / (2 * ny));
	jy = jy.sin().square();
	jy = jy * (-4 / (h_y * h_y));

	RowVectorXd x = jx;
	VectorXd    y = jy;
	denom         = VectorXd::Ones(ny) * x + y * RowVectorXd::Ones(nx);
}
void Domain::solve()
{
	// cout << "North:\n" << boundary_north << "\n\n";
	// cout << "South:\n" << boundary_south << "\n\n";
	// cout << "East:\n" << boundary_east << "\n\n";
	// cout << "West:\n" << boundary_west << "\n\n";
	MyArray grid_copy = grid;
	grid_copy.row(0) += -2 / ((h_y) * (h_y)) * boundary_north;
	grid_copy.row(ny - 1) += -2 / ((h_y) * (h_y)) * boundary_south;
	grid_copy.col(0) += -2 / ((h_x) * (h_x)) * boundary_west;
	grid_copy.col(nx - 1) += -2 / ((h_x) * (h_x)) * boundary_east;

	MyArray tmp = MyArray(ny, nx);

	fftw_plan plan1
	= fftw_plan_r2r_2d(ny, nx, &grid_copy(0), &tmp(0), FFTW_RODFT10, FFTW_RODFT10, FFTW_MEASURE);

	fftw_plan plan2
	= fftw_plan_r2r_2d(ny, nx, &tmp(0), &u(0), FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);

	fftw_execute(plan1);
	tmp /= denom;
	fftw_execute(plan2);
	u /= 4 * nx * ny;
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
}
