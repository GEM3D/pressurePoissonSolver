#include <Eigen/Dense>
typedef Eigen::Array<double, -1, -1> MyArray;
class Domain
{
	public:
	MyArray boundary_north;
	MyArray boundary_south;
	MyArray boundary_east;
	MyArray boundary_west;
	MyArray grid;
	MyArray denom;
	MyArray u;
	double  h_x;
	double  h_y;
	int     nx;
	int     ny;
	Domain(MyArray grid, double h_x, double h_y);
	Domain() {}
	void solve();
};
