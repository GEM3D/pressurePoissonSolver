#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/SparseCholesky>
class Domain
{
	public:
	Eigen::RowVectorXd                                  boundary_north;
	Eigen::RowVectorXd                                  boundary_south;
	Eigen::VectorXd                                     boundary_east;
	Eigen::VectorXd                                     boundary_west;
	Eigen::MatrixXd                                     grid;
	Eigen::MatrixXd                                     u;
	double                                              h_x;
	double                                              h_y;
	Eigen::SparseMatrix<double>                         A;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> *solver;
	Domain(Eigen::SparseMatrix<double> A, Eigen::MatrixXd grid,double h_x,double h_y);
	void solve();
};
