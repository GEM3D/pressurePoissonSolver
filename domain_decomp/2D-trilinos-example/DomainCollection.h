#ifndef DOMAINCOLLECTION_H
#define DOMAINCOLLECTION_H
#include "Domain.h"
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
#include "PatchSolver.h"
#include "Interpolator.h"
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
	Teuchos::RCP<PatchSolver> solver;
	bool                      neumann = false;
	bool                      zero_u  = false;
	/**
	 * @brief A map of domain ids to the domain objects.
	 */
	/**
	 * @brief The MPI communicator used.
	 */
	Teuchos::RCP<const Teuchos::Comm<int>> comm;
	int                                    num_cols = 0;

	/**
	 * @brief Number of global domains
	 */
	int num_global_domains;

	public:
    Teuchos::RCP<Interpolator> interpolator;
	void setPatchSolver(Teuchos::RCP<PatchSolver> psolver);
	/**
	 * @brief The number of cells in a single dimension.
	 *
	 * i.e. n means nxn cells in each domain.
	 */
	int n;

	std::map<int, Teuchos::RCP<Domain>> domains;
	bool   amr    = false;
	double f_mean = 0;
	/**
	 * @brief Tpetra map that assures each domains has access to it's interface values
	 */
	/**
	 * @brief Tpetra map that assures each domains has access to it's interface information values
	 */
	Teuchos::RCP<map_type> collection_iface_map;
	Teuchos::RCP<map_type> collection_row_iface_map;
	/**
	 * @brief Tpetra map use for interface information.
	 */
	Teuchos::RCP<map_type>    iface_map;
	Teuchos::RCP<map_type>    row_iface_map;
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
	DomainCollection(DomainSignatureCollection dsc, int n,
	                 Teuchos::RCP<const Teuchos::Comm<int>> comm);

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
	void solveWithInterface(const vector_type &f, vector_type &u, const vector_type &gamma,
	                        vector_type &diff);

	/**
	 * @brief get the difference in flux on refined boundaries
	 */
	Teuchos::RCP<vector_type> getInterfaceCoords();

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
	double integrateBoundaryFlux();
	double area();
	double integrateAU();
	void   setZeroU() { zero_u = true; }
	/**
	 * @brief Form the Schur complement matrix using an RBMatrix
	 *
	 * @param map the map used in the matrix
	 *
	 * @return the formed matrix
	 */

	typedef std::function<void(int, int, Teuchos::RCP<std::valarray<double>>, bool, bool)> inserter;
	void assembleMatrix(inserter insertBlock, int n = -1);
	void formCRSMatrix(Teuchos::RCP<map_type> map, Teuchos::RCP<matrix_type> &A, int n = -1);

	void swapResidSol()
	{
		for (auto &p : domains) {
			p.second->swapResidSol();
		}
	}
	void sumResidIntoSol()
	{
		for (auto &p : domains) {
			p.second->sumResidIntoSol();
		}
	}
	void outputClaw();
#ifdef HAVE_VTK
	void outputVTK();
#endif
	int getGlobalNumCells() { return num_global_domains * n * n; }
};
#endif
