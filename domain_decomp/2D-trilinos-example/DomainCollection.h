#ifndef DOMAINCOLLECTION_H
#define DOMAINCOLLECTION_H
#include "Domain.h"
#include "Iface.h"
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
#include "RBMatrix.h"
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
	/**
	 * @brief A map of domain ids to the domain objects.
	 */
	std::map<int, Teuchos::RCP<Domain>> domains;
	/**
	 * @brief The MPI communicator used.
	 */
	Teuchos::RCP<const Teuchos::Comm<int>> comm;
	/**
	 * @brief The number of cells in a single dimension.
	 *
	 * i.e. n means nxn cells in each domain.
	 */
	int n;

	int num_cols = 0;

	/**
	 * @brief Number of global domains
	 */
	int num_global_domains;
	/**
     * TODO update this
	 * @brief a vector of integers that stores information for each inteface.
	 *
	 * For the current format each interface requires 15 integers in order to store information.
	 *
	 * Each index represents:
	 *   0. The interface's global index in the interface value array.
	 *   1. The interfaces type
	 *   2. The interfaces western neighbor on it's right
	 *   3. The interfaces western neighbor on it's right type
	 *   4. The interfaces northern neighbor on it's right
	 *   5. The interfaces northern neighbor on it's right type
	 *   6. The interfaces eastern neighbor on it's right
	 *   7. The interfaces eastern neighbor on it's right type
	 *   8. The interfaces western neighbor on it's left
	 *   9. The interfaces western neighbor on it's left type
	 *   10. The interfaces northern neighbor on it's left
	 *   11. The interfaces northern neighbor on it's left type
	 *   12. The interfaces eastern neighbor on it's left
	 *   13. The interfaces eastern neighbor on it's left type
	 *   14. The axis that the interface resides on
	 */
	Teuchos::RCP<int_vector_type> iface_info;

    std::set<Iface> ifaces;

	public:
	bool amr = false;
	double f_mean=0;
	/**
	 * @brief Tpetra map that assures each domains has access to it's interface values
	 */
	Teuchos::RCP<map_type> collection_map;
	/**
	 * @brief Tpetra map that assures each domains has access to it's interface information values
	 */
	Teuchos::RCP<map_type> collection_iface_map;
	/**
	 * @brief Tpetra map for interface values used in Schur complement matrix
	 */
	Teuchos::RCP<map_type> matrix_map;
	/**
	 * @brief Tpetra map use for interface information.
	 */
	Teuchos::RCP<map_type> iface_map;

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
	void solveWithInterface(const vector_type &gamma, vector_type &diff);
	void solveWithInterface(const vector_type &gamma);

    /**
     * @brief get the difference in flux on refined boundaries
     */
	void getFluxDiff(vector_type &diff);

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
	double uSum();

	/**
	 * @return sum of exact solution
	 */
	double exactSum();

	/**
	 * @return the residual
	 */
	double residual();

	double integrateF();
	double integrateAU();

	/**
	 * @brief Form the Schur complement matrix using an RBMatrix
	 *
	 * @param map the map used in the matrix
	 *
	 * @return the formed matrix
	 */
	Teuchos::RCP<RBMatrix> formRBMatrix(Teuchos::RCP<map_type> map, int delete_row = -1);

	void swapResidSol(){
		for (auto &p : domains) {
			p.second->swapResidSol();
		}
	}
	void sumResidIntoSol(){
		for (auto &p : domains) {
			p.second->sumResidIntoSol();
		}
    }
	void outputClaw();
	int  getGlobalNumCells() { return num_global_domains * n * n; }
};
#endif
