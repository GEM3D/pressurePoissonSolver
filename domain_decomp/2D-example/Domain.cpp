#include "Domain.h"
using namespace Eigen;
Domain::Domain(SparseMatrix<double> A,MatrixXd grid,double h_x, double h_y){
    this->A = A;
    this->grid = grid;
    this->u = MatrixXd(grid.rows(),grid.cols());
    solver = new SimplicialLDLT<SparseMatrix<double>>();
    solver->compute(A);
    this->h_x = h_x;
    this->h_y = h_y;
    boundary_north = RowVectorXd(grid.cols());
    boundary_south = RowVectorXd(grid.cols());
}
void Domain::solve(){
    MatrixXd grid_copy = grid;
    grid_copy.row(0) += -2/((h_y)*(h_y))*boundary_north;
    grid_copy.row(grid.rows()-1) += -2/((h_y)*(h_y))*boundary_south;
    grid_copy.col(0) += -2/((h_x)*(h_x))*boundary_east;
    grid_copy.col(grid.cols()-1) += -2/((h_x)*(h_x))*boundary_west;
    Map<VectorXd> f(grid_copy.data(),grid_copy.size());
    VectorXd sol = solver->solve(f);
    u = Map<MatrixXd>(sol.data(),grid.rows(),grid.cols());
}
