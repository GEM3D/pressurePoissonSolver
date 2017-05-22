#ifndef DOMAINCOLLECTION_H
#define DOMAINCOLLECTION_H
#include "Domain.h"
#include "DomainSignatureCollection.h"
#include <HYPRE_sstruct_ls.h>
#include <functional>
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
	int                   n;
	int                   num_global_domains;
	bool                  neumann    = false;
	bool                  use_parcsr = false;
	HYPRE_SStructGrid     grid;
	HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
	HYPRE_SStructGraph    graph;
	HYPRE_SStructStencil  stencil_5pt;

	public:
	HYPRE_SStructMatrix   A;
	HYPRE_SStructVector   b;
	HYPRE_SStructVector   x;
	HYPRE_ParCSRMatrix    par_A;
	HYPRE_ParVector       par_b;
	HYPRE_ParVector       par_x;
	std::map<int, Domain> domains;
	bool amr = false;
	double f_mean=0;
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
	DomainCollection(DomainSignatureCollection dsc, int n);

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
	void saveResult();
	void initVectors();
	void formMatrix();

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

	void swapResidSol(){
		for (auto &p : domains) {
			p.second.swapResidSol();
		}
	}
	void sumResidIntoSol(){
		for (auto &p : domains) {
			p.second.sumResidIntoSol();
		}
    }
	void outputClaw();
	int  getGlobalNumCells() { return num_global_domains * n * n; }
	void setParCSR() { use_parcsr = true; }
};
#endif
