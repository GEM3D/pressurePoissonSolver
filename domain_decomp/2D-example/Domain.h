#ifndef DOMAIN_H
#define DOMAIN_H
#include <Eigen/Dense>
#include <fftw3.h>
typedef Eigen::Array<double, -1, -1> MyArray;

class Domain;
class Interface
{
	public:
	enum class axis { x, y };
	int             start_index;
	int             size;
	axis            dir;
	Domain *        left  = nullptr;
	Domain *        right = nullptr;
	Eigen::VectorXd left_zero;
	Eigen::VectorXd right_zero;
	MyArray         gamma;
	Interface() {}
	Interface(int start_index, int size, axis dir)
	{
		this->start_index = start_index;
		this->size        = size;
		this->dir         = dir;
		if (dir == axis::x) {
			gamma = Eigen::RowVectorXd(size);
		} else {
			gamma = Eigen::VectorXd(size);
		}
	}
	MyArray getDiff();
};
class Domain
{
	public:
	//	MyArray boundary_north;
	//	MyArray boundary_south;
	//	MyArray boundary_east;
	//	MyArray boundary_west;
	MyArray    grid;
	MyArray    grid_copy;
	MyArray    tmp;
	MyArray    denom;
	MyArray    u;
	Interface *north        = nullptr;
	Interface *east         = nullptr;
	Interface *south        = nullptr;
	Interface *west         = nullptr;
	Domain *   north_domain = nullptr;
	Domain *   east_domain  = nullptr;
	Domain *   south_domain = nullptr;
	Domain *   west_domain  = nullptr;
	fftw_plan  plan1;
	fftw_plan  plan2;
	double     h_x;
	double     h_y;
	int        nx;
	int        ny;
	Domain(MyArray grid, double h_x, double h_y);
	Domain() {}
	void solve();
};
#endif
