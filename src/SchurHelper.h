/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively 
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef SCHURHELPER_H
#define SCHURHELPER_H
#include "DomainCollection.h"
#include "Iface.h"
#include "Interpolator.h"
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
template <size_t D> class SchurHelper
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
	std::shared_ptr<Interpolator<D>> interpolator;

	/**
	 * @brief The patch operator
	 */
	std::shared_ptr<PatchOperator<D>> op;

	/**
	 * @brief The patch solver
	 */
	std::shared_ptr<PatchSolver<D>> solver;

	std::deque<SchurDomain<D>> domains;
	std::map<int, IfaceSet<D>> ifaces;

	std::vector<int> iface_dist_map_vec;
	std::vector<int> iface_map_vec;
	std::vector<int> iface_off_proc_map_vec;
	void             indexIfacesLocal();
	void             indexDomainIfacesLocal();
	void             indexIfacesGlobal();

	int num_global_ifaces = 0;

	public:
	SchurHelper() = default;
	/**
	 * @brief Create a SchurHelper from a given DomainCollection
	 *
	 * @param dc the DomainCollection
	 * @param comm the teuchos communicator
	 */
	SchurHelper(DomainCollection<D> dc, std::shared_ptr<PatchSolver<D>> solver,
	            std::shared_ptr<PatchOperator<D>> op,
	            std::shared_ptr<Interpolator<D>>  interpolator);

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

	PW_explicit<Vec> getNewSchurVec()
	{
		PW<Vec> u;
		VecCreateMPI(MPI_COMM_WORLD, iface_map_vec.size() * std::pow(n, D - 1), PETSC_DETERMINE,
		             &u);
		return u;
	}
	PW_explicit<Vec> getNewSchurDistVec()
	{
		PW<Vec> u;
		VecCreateSeq(PETSC_COMM_SELF, iface_dist_map_vec.size() * std::pow(n, D - 1), &u);
		return u;
	}

	int getSchurVecLocalSize()
	{
		return iface_map_vec.size() * std::pow(n, D - 1);
	}
	int getSchurVecGlobalSize()
	{
		return num_global_ifaces * std::pow(n, D - 1);
	}
	// getters
	std::shared_ptr<Interpolator<D>> getInterpolator()
	{
		return interpolator;
	}
	std::shared_ptr<PatchOperator<D>> getOp()
	{
		return op;
	}
	std::shared_ptr<PatchSolver<D>> getSolver()
	{
		return solver;
	}
	const std::map<int, IfaceSet<D>> getIfaces() const
	{
		return ifaces;
	}
	int getN() const
	{
		return n;
	}
};
template <size_t D>
inline SchurHelper<D>::SchurHelper(DomainCollection<D> dc, std::shared_ptr<PatchSolver<D>> solver,
                                   std::shared_ptr<PatchOperator<D>> op,
                                   std::shared_ptr<Interpolator<D>>  interpolator)
{
	this->n = dc.getN();
	for (auto &p : dc.domains) {
		domains.push_back(*p.second);
	}
	std::map<int, std::pair<int, IfaceSet<D>>> off_proc_ifaces;
	for (SchurDomain<D> &sd : domains) {
		sd.enumerateIfaces(ifaces, off_proc_ifaces);
		solver->addDomain(sd);
	}
	{
        using namespace std;
	    // send info
	    deque<char *>       buffers;
	    deque<char *>       recv_buffers;
	    vector<MPI_Request> requests;
	    vector<MPI_Request> send_requests;
	    for (auto &p : off_proc_ifaces) {
	        int       dest   = p.second.first;
	        IfaceSet<D> &iface  = p.second.second;
	        int       size   = iface.serialize(nullptr);
	        char *    buffer = new char[size];
	        buffers.push_back(buffer);
	        iface.serialize(buffer);
	        MPI_Request request;
	        MPI_Isend(buffer, size, MPI_CHAR, dest, 0, MPI_COMM_WORLD, &request);
	        send_requests.push_back(request);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
	    int        is_message;
	    MPI_Status status;
	    MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &is_message, &status);
	    // recv info
	    while (is_message) {
	        int size;
	        MPI_Get_count(&status, MPI_CHAR, &size);
	        char *buffer = new char[size];
	        recv_buffers.push_back(buffer);

	        MPI_Request request;
	        MPI_Irecv(buffer, size, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
	                  &request);
	        MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &is_message, &status);
	        requests.push_back(request);
	    }
	    // wait for all
	    MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Startall(send_requests.size(), &send_requests[0]);
	    MPI_Startall(requests.size(), &requests[0]);
	    MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
	    MPI_Barrier(MPI_COMM_WORLD);
	    // delete send buffers
	    for (char *buffer : buffers) {
	        delete[] buffer;
	    }
	    // process received objects
	    for (char *buffer : recv_buffers) {
	        IfaceSet<D> ifs;
	        ifs.deserialize(buffer);
	        ifaces[ifs.id].insert(ifs);
	        delete[] buffer;
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
	}
	indexDomainIfacesLocal();
	indexIfacesLocal();
	this->solver       = solver;
	this->op           = op;
	this->interpolator = interpolator;
	local_gamma        = getNewSchurDistVec();
	local_interp       = getNewSchurDistVec();
	gamma              = getNewSchurVec();
	PW<IS> dist_is;
	ISCreateBlock(MPI_COMM_SELF, std::pow(n, D - 1), iface_dist_map_vec.size(),
	              &iface_dist_map_vec[0], PETSC_COPY_VALUES, &dist_is);
	VecScatterCreate(gamma, dist_is, local_gamma, nullptr, &scatter);

	int num_ifaces = ifaces.size();
	MPI_Allreduce(&num_ifaces, &num_global_ifaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}
template <size_t D>
inline void SchurHelper<D>::solveWithInterface(const Vec f, Vec u, const Vec gamma, Vec diff)
{
	// initilize our local variables
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	VecScale(local_interp, 0);
	// solve over domains on this proc
	solver->domainSolve(domains, f, u, local_gamma);
	for (SchurDomain<D> &sd : domains) {
		interpolator->interpolate(sd, u, local_interp);
	}

	// export diff vector
	VecScale(diff, 0);
	VecScatterBegin(scatter, local_interp, diff, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, diff, ADD_VALUES, SCATTER_REVERSE);
	VecAXPBY(diff, 1.0, -1.0, gamma);
}
template <size_t D>
inline void SchurHelper<D>::solveAndInterpolateWithInterface(const Vec f, Vec u, const Vec gamma,
                                                             Vec interp)
{
	// initilize our local variables
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	VecScale(local_interp, 0);
	// solve over domains on this proc
	solver->domainSolve(domains, f, u, local_gamma);
	for (SchurDomain<D> &sd : domains) {
		interpolator->interpolate(sd, u, local_interp);
	}

	// export diff vector
	VecScale(interp, 0);
	VecScatterBegin(scatter, local_interp, interp, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, interp, ADD_VALUES, SCATTER_REVERSE);
}
template <size_t D> inline void SchurHelper<D>::solveWithSolution(const Vec f, Vec u)
{
	// initilize our local variables
	VecScale(local_gamma, 0);
	VecScale(local_interp, 0);
	for (SchurDomain<D> &sd : domains) {
		interpolator->interpolate(sd, u, local_interp);
	}
	VecScale(gamma, 0);
	VecScatterBegin(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	// solve over domains on this proc
	solver->domainSolve(domains, f, u, local_gamma);
}
template <size_t D>
inline void SchurHelper<D>::interpolateToInterface(const Vec f, Vec u, Vec gamma)
{
	// initilize our local variables
	VecScale(local_interp, 0);
	for (SchurDomain<D> &sd : domains) {
		interpolator->interpolate(sd, u, local_interp);
	}
	VecScale(gamma, 0);
	VecScatterBegin(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
}
template <size_t D>
inline void SchurHelper<D>::applyWithInterface(const Vec u, const Vec gamma, Vec f)
{
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	for (SchurDomain<D> &sd : domains) {
		op->apply(sd, u, local_gamma, f);
	}
}
template <size_t D> inline void SchurHelper<D>::apply(const Vec u, Vec f)
{
	VecScale(local_interp, 0);
	for (SchurDomain<D> &sd : domains) {
		interpolator->interpolate(sd, u, local_interp);
	}
	VecScale(gamma, 0);
	VecScatterBegin(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	for (SchurDomain<D> &sd : domains) {
		op->apply(sd, u, local_gamma, f);
	}
}
template <size_t D> inline void SchurHelper<D>::indexDomainIfacesLocal()
{
	using namespace std;
	vector<int>   map_vec;
	map<int, int> rev_map;
	if (!domains.empty()) {
		int curr_i = 0;
		for (SchurDomain<D> &sd : domains) {
			for (int id : sd.getIds()) {
				if (rev_map.count(id) == 0) {
					rev_map[id] = curr_i;
					map_vec.push_back(id);
					curr_i++;
				}
			}
		}
		for (SchurDomain<D> &sd : domains) {
			sd.setLocalIndexes(rev_map);
		}
	}
	iface_dist_map_vec = map_vec;
}
template <size_t D> inline void SchurHelper<D>::indexIfacesLocal()
{
	using namespace std;
	int           curr_i = 0;
	vector<int>   map_vec;
	vector<int>   off_proc_map_vec;
	vector<int>   off_proc_map_vec_send;
	map<int, int> rev_map;
	if (!ifaces.empty()) {
		set<int> todo;
		for (auto &p : ifaces) {
			todo.insert(p.first);
		}
		set<int> enqueued;
		while (!todo.empty()) {
			deque<int> queue;
			queue.push_back(*todo.begin());
			enqueued.insert(*todo.begin());
			while (!queue.empty()) {
				int i = queue.front();
				todo.erase(i);
				queue.pop_front();
				map_vec.push_back(i);
				IfaceSet<D> &ifs = ifaces.at(i);
				rev_map[i]       = curr_i;
				curr_i++;
				for (int nbr : ifs.getNbrs()) {
					if (!enqueued.count(nbr)) {
						enqueued.insert(nbr);
						if (ifaces.count(nbr)) {
							queue.push_back(nbr);
						} else {
							off_proc_map_vec.push_back(nbr);
						}
					}
				}
			}
		}
	}

	// map off proc
	for (int i : off_proc_map_vec) {
		rev_map[i] = curr_i;
		curr_i++;
	}
	for (auto &p : ifaces) {
		p.second.setLocalIndexes(rev_map);
	}
	iface_map_vec          = map_vec;
	iface_off_proc_map_vec = off_proc_map_vec;
	indexIfacesGlobal();
}
template <size_t D> inline void SchurHelper<D>::indexIfacesGlobal()
{
	using namespace std;
	// global indices are going to be sequentially increasing with rank
	int local_size = ifaces.size();
	int start_i;
	MPI_Scan(&local_size, &start_i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	start_i -= local_size;
	vector<int> new_global(local_size);
	iota(new_global.begin(), new_global.end(), start_i);

	// create map for gids
	PW<AO> ao;
	AOCreateMapping(MPI_COMM_WORLD, local_size, &iface_map_vec[0], &new_global[0], &ao);

	// get indices for schur matrix
	{
		// get global indices that we want to recieve for dest vector
		vector<int> inds = iface_map_vec;
		for (int i : iface_off_proc_map_vec) {
			inds.push_back(i);
		}

		// get new global indices
		AOApplicationToPetsc(ao, inds.size(), &inds[0]);
		map<int, int> rev_map;
		for (size_t i = 0; i < inds.size(); i++) {
			rev_map[i] = inds[i];
		}

		// set new global indices in iface objects
		for (auto &p : ifaces) {
			p.second.setGlobalIndexes(rev_map);
		}
		for (size_t i = 0; i < iface_map_vec.size(); i++) {
			iface_map_vec[i] = inds[i];
		}
		for (size_t i = 0; i < iface_off_proc_map_vec.size(); i++) {
			iface_off_proc_map_vec[i] = inds[iface_map_vec.size() + i];
		}
	}
	// get indices for local ifaces
	{
		// get global indices that we want to recieve for dest vector
		vector<int> inds = iface_dist_map_vec;

		// get new global indices
		AOApplicationToPetsc(ao, inds.size(), &inds[0]);
		map<int, int> rev_map;
		for (size_t i = 0; i < inds.size(); i++) {
			rev_map[i] = inds[i];
		}

		// set new global indices in domain objects
		for (SchurDomain<D> &sd : domains) {
			sd.setGlobalIndexes(rev_map);
		}
		for (size_t i = 0; i < iface_dist_map_vec.size(); i++) {
			iface_dist_map_vec[i] = inds[i];
		}
	}
}
#endif
