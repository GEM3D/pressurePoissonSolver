#ifndef DOMAINCOLLECTION_H
#define DOMAINCOLLECTION_H
#include "DomainSignatureCollection.h"
#include <functional>
class Domain;
struct TrilinosObjects;
struct AmgxCrs{
    int n_global;
    int n;
    int nnz;
    int block_dimx=1;
    int block_dimy=1;
    std::vector<int> row_ptrs;
    std::vector<int64_t>  col_indicies_global;
    std::vector<double> data;
  double * diag_data=nullptr;
};
/**
 * @brief This class represents a collection of domains that a single processor owns.
 *
 * The purposes of this class:
 *   - Provide a member function for solving with a given interface vector.
 *   - Handle the initialization of the domains.
 *   - Provide member functions for calculating error, residual, etc.
 *   - Provide member functions that generate the Schur complement matrix.
 */
class DomainCollection
{
	private:
	bool neumann = false;
	bool zero_u  = false;
	int                                    num_cols = 0;

	/**
	 * @brief Number of global domains
	 */
	int num_global_domains;
    TrilinosObjects * trilObjs;

	public:
	/**
	 * @brief The number of cells in a single dimension.
	 *
	 * i.e. n means nxn cells in each domain.
	 */
	int n;

	std::map<int, Domain*> domains;
	bool   amr    = false;
	double f_mean = 0;
	/**
	 * @brief Tpetra map that assures each domains has access to it's interface values
	 */

	DomainSignatureCollection dsc;

	/**
	 * @brief Create a DomainCollection from a given DomainSignatureCollection
	 *
	 * @param dsc the DomainSignatureCollection
	 * @param n number of cells in each direction for each domain
	 * @param h_x the x spacing
	 * @param h_y the y spacing
	 * @param comm the teuchos communicator
	 */
	DomainCollection(DomainSignatureCollection dsc, int n
	                 );

	/**
	 * @brief Initialize domains using Dirichlet boundary conditions.
	 *
	 * @param ffun the function we are solving for
	 * @param gfun the exact solution
	 */
	void initDirichlet(std::function<double(double, double)> ffun,
	                   std::function<double(double, double)> gfun);
	/**
	 * @brief Initialize domains using Neumann boundary conditions
	 *
	 * @param ffun the function we are solving for
	 * @param efun the exact solution
	 * @param nfunx the d/dx derivative
	 * @param nfuny the d/dy derivative
	 */
	void initNeumann(std::function<double(double, double)> ffun,
	                 std::function<double(double, double)> efun,
	                 std::function<double(double, double)> nfunx,
	                 std::function<double(double, double)> nfuny, bool amr);

	/**
	 * @brief output solution in MatrixMarket format
	 *
	 * @param out the stream to output to
	 */
	void outputSolution(std::ostream &out);
	void outputSolutionRefined(std::ostream &out);

	/**
	 * @brief output residual in MatrixMarket format
	 *
	 * @param out the stream to output to
	 */
	void outputResidual(std::ostream &out);
	void outputResidualRefined(std::ostream &out);

	/**
	 * @brief output error in MatrixMarket format
	 *
	 * @param out the stream to output to
	 */
	void outputError(std::ostream &out);
	void outputErrorRefined(std::ostream &out);

	/**
	 * @brief Solve with a given set of interface values
	 *
	 * @param gamma the interface values to use
	 * @param diff the resulting difference
	 */
	void solveWithInterface(double* gamma, double * diff);

	/**
	 * @brief Generate Tpetra maps
	 */
	void generateMaps();

	/**
	 * @brief Distribute interface information
	 */
	void distributeIfaceInfo();

	/**
	 * @return norm of difference of exact and computed solution
	 */
	double diffNorm();

	/**
	 * @brief norm function used in Neumann case
	 *
	 * @param uavg average of computed solution
	 * @param eavg average of exact solution
	 *
	 * @return norm of difference of exact and computed solution
	 */
	double diffNorm(double uavg, double eavg);

	/**
	 * @return norm of exact solution
	 */
	double exactNorm();

	/**
	 * @return norm of right hand side
	 */
	double fNorm();

	/**
	 * @brief norm function used in Neumann case
	 *
	 * @param eavg average of exact solution
	 *
	 * @return norm of exact solution
	 */
	double exactNorm(double eavg);

	/**
	 * @return sum of computed solution
	 */
	double integrateU();

	/**
	 * @return sum of exact solution
	 */
	double integrateExact();

	/**
	 * @return the residual
	 */
	double residual();

	double integrateF();
	void zeroF(double fdiff);
	double integrateBoundaryFlux();
	double area();
	double integrateAU();
	void   setZeroU() { zero_u = true; }
	AmgxCrs formCRSMatrix();

	void swapResidSol();
	void sumResidIntoSol();
	void outputClaw();
#ifdef HAVE_VTK
	void outputVTK();
#endif
	int  getGlobalNumCells() { return num_global_domains * n * n; }
};
#endif
