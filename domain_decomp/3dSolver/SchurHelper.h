#ifndef SCHURHELPER_H
#define SCHURHELPER_H
#include "DomainCollection.h"
#include "Iface.h"
#include "Interpolator.h"
#include "PBMatrix.h"
#include "PatchOperator.h"
#include "PatchSolvers/PatchSolver.h"
#include "SchurDomain.h"
#include <deque>
#include <memory>
#include <petscmat.h>
#include <petscpc.h>
#include <valarray>
/**
 * @brief This class represents a collection of domains that a single processor owns.
 *
 * The purposes of this class:
 *   - Provide a member function for solving with a given interface vector.
 *   - Handle the initialization of the domains.
 *   - Provide member functions for calculating error, residual, etc.
 *   - Provide member functions that generate the Schur complement matrix.
 */
struct Block;
class SchurHelper
{
	private:
	int n;

	PW<Vec>        local_gamma;
	PW<Vec>        gamma;
	PW<Vec>        local_interp;
	PW<VecScatter> scatter;

	/**
	 * @brief Interpolates to interface values
	 */
	std::shared_ptr<Interpolator> interpolator;

	/**
	 * @brief The patch operator
	 */
	std::shared_ptr<PatchOperator> op;

	/**
	 * @brief The patch solver
	 */
	std::shared_ptr<PatchSolver> solver;

	typedef std::function<void(Block *, std::shared_ptr<std::valarray<double>>)> inserter;
	void assembleMatrix(inserter insertBlock);

	std::deque<SchurDomain<3>> domains;
	std::map<int, IfaceSet>    ifaces;

	std::vector<int> iface_dist_map_vec;
	std::vector<int> iface_map_vec;
	std::vector<int> iface_off_proc_map_vec;
	void             indexIfacesLocal();
	void             indexDomainIfacesLocal();
	void             indexIfacesGlobal();
	void             zoltanBalance();

	int num_global_ifaces = 0;

	public:
	SchurHelper() = default;
	/**
	 * @brief Create a SchurHelper from a given DomainCollection
	 *
	 * @param dc the DomainCollection
	 * @param comm the teuchos communicator
	 */
	SchurHelper(DomainCollection dc, std::shared_ptr<PatchSolver> solver,
	            std::shared_ptr<PatchOperator> op, std::shared_ptr<Interpolator> interpolator);

	/**
	 * @brief Solve with a given set of interface values
	 *
	 * @param f the rhs vector
	 * @param u the vector to put solution in
	 * @param gamma the interface values to use
	 * @param diff the resulting difference
	 */
	void solveWithInterface(const Vec f, Vec u, const Vec gamma, Vec diff);
	void solveAndInterpolateWithInterface(const Vec f, Vec u, const Vec gamma, Vec interp);
	void solveWithSolution(const Vec f, Vec u);
	void interpolateToInterface(const Vec f, Vec u, Vec gamma);

	/**
	 * @brief Apply patch operator with a given set of interface values
	 *
	 * @param u the solution vector to use
	 * @param gamma the interface values to use
	 * @param f the resulting rhs vector
	 */
	void applyWithInterface(const Vec u, const Vec gamma, Vec f);
	void apply(const Vec u, Vec f);

	/**
	 * @brief Form the Schur complement matrix
	 *
	 * @return the formed matrix
	 */
	PW_explicit<Mat> formCRSMatrix();
	PBMatrix *       formPBMatrix();
	void             getPBDiagInv(PC p);
	PW_explicit<Mat> getPBMatrix();
	PW_explicit<Mat> getPBDiagInv();
	PW_explicit<Vec> getNewSchurVec();
	PW_explicit<Vec> getNewSchurDistVec();

	int getSchurVecLocalSize()
	{
		return iface_map_vec.size() * n * n;
	}
	int getSchurVecGlobalSize()
	{
		return num_global_ifaces * n * n;
	}
	// getters
	std::shared_ptr<Interpolator> getInterpolator()
	{
		return interpolator;
	}
	std::shared_ptr<PatchOperator> getOp()
	{
		return op;
	}
	std::shared_ptr<PatchSolver> getSolver()
	{
		return solver;
	}
};
#endif
