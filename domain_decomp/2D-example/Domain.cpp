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
	jx         = jx * (M_PI / (2 * nx));
	jx         = jx.sin().square();
	jx         = jx * (-4 / (h_x * h_x));

	MyArray jy = VectorXd::LinSpaced(ny, 1, ny);
	jy         = jy * (M_PI / (2 * ny));
	jy         = jy.sin().square();
	jy         = jy * (-4 / (h_y * h_y));

	RowVectorXd x = jx;
	VectorXd    y = jy;
	denom         = VectorXd::Ones(ny) * x + y * RowVectorXd::Ones(nx);
	grid_copy     = grid;
	tmp           = MyArray(ny, nx);
	plan1
	= fftw_plan_r2r_2d(nx, ny, &grid_copy(0), &tmp(0), FFTW_RODFT10, FFTW_RODFT10, FFTW_MEASURE);

	plan2 = fftw_plan_r2r_2d(nx, ny, &tmp(0), &u(0), FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);
}
void Domain::solve()
{
	// cout << "North:\n" << boundary_north << "\n\n";
	// cout << "South:\n" << boundary_south << "\n\n";
	// cout << "East:\n" << boundary_east << "\n\n";
	// cout << "West:\n" << boundary_west << "\n\n";
	grid_copy = grid;
	if (north) {
		grid_copy.row(ny - 1) += -2 / ((h_y) * (h_y)) * north->gamma;
	}
	if (east) {
		grid_copy.col(nx - 1) += -2 / ((h_x) * (h_x)) * east->gamma;
	}
	if (south) {
		grid_copy.row(0) += -2 / ((h_y) * (h_y)) * south->gamma;
	}
	if (west) {
		grid_copy.col(0) += -2 / ((h_x) * (h_x)) * west->gamma;
	}

	fftw_plan plan1
	= fftw_plan_r2r_2d(nx, ny, &grid_copy(0), &tmp(0), FFTW_RODFT10, FFTW_RODFT10, FFTW_MEASURE);

	fftw_plan plan2
	= fftw_plan_r2r_2d(nx, ny, &tmp(0), &u(0), FFTW_RODFT01, FFTW_RODFT01, FFTW_MEASURE);

	fftw_execute(plan1);
	tmp /= denom;
	fftw_execute(plan2);
	u /= 4 * nx * ny;
}

MyArray Interface::getDiff()
{
	if (dir == axis::x) {
		return (left->u.row(left->u.rows() - 1) + right->u.row(0) - 2 * gamma).transpose();
	} else {
		return left->u.col(left->u.cols() - 1) + right->u.col(0) - 2 * gamma;
	}
}

MyArray Interface::getZeroDiffFromLeft()
{
	if (dir == axis::x) {
		return (left->u.row(left->u.rows() - 1) + right_zero - 2 * gamma).transpose();
	} else {
		return left->u.col(left->u.cols() - 1) + right_zero - 2 * gamma;
	}
}

MyArray Interface::getZeroDiffFromRight()
{
	if (dir == axis::x) {
		return (left_zero + right->u.row(0) - 2 * gamma).transpose();
	} else {
		return left_zero + right->u.col(0) - 2 * gamma;
	}
}

void Interface::setZero()
{
	if (dir == axis::x) {
		left_zero  = left->u.row(left->u.rows() - 1);
		right_zero = right->u.row(0);
	} else {
		left_zero  = left->u.col(left->u.cols() - 1);
		right_zero = right->u.col(0);
	}
}
