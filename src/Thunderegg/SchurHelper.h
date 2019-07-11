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

#ifndef THUNDEREGG_SCHURHELPER_H
#define THUNDEREGG_SCHURHELPER_H
#include <Thunderegg/Domain.h>
#include <Thunderegg/Iface.h>
#include <Thunderegg/IfaceInterp.h>
#include <Thunderegg/PatchOperator.h>
#include <Thunderegg/PatchSolvers/PatchSolver.h>
#include <Thunderegg/PetscVector.h>
#include <Thunderegg/SchurInfo.h>
#include <deque>
#include <memory>
#include <petscmat.h>
#include <petscpc.h>
#include <valarray>
/**
 * @brief This class represents a collection of sinfo_vector that a single processor owns.
 *
 * The purposes of this class:
 *   - Provide a member function for solving with a given interface vector.
 */
template <size_t D> class SchurHelper
{
	private:
	std::shared_ptr<Domain<D>> domain;

	std::shared_ptr<PetscVector<D - 1>> local_gamma;
	std::shared_ptr<PetscVector<D - 1>> gamma;
	PW<VecScatter>                      scatter;

	/**
	 * @brief Interpolates to interface values
	 */
	std::shared_ptr<IfaceInterp<D>> interpolator;
	/**
	 * @brief The patch operator
	 */
	std::shared_ptr<PatchOperator<D>> op;
	/**
	 * @brief The patch solver
	 */
	std::shared_ptr<PatchSolver<D>> solver;
	/**
	 * @brief Vector of SchurInfo pointers where index in the vector corresponds to the patch's
	 * local index
	 */
	std::vector<SchurInfo<D>>  sinfo_vector;
	std::map<int, IfaceSet<D>> ifaces;

	std::vector<int> iface_dist_map_vec;
	std::vector<int> iface_map_vec;
	std::vector<int> iface_off_proc_map_vec;
	void             indexIfacesLocal();
	void             indexDomainIfacesLocal();
	void             indexIfacesGlobal();

	int                    num_global_ifaces = 0;
	int                    iface_stride;
	std::array<int, D - 1> lengths;

	public:
	SchurHelper() = default;
	/**
	 * @brief Create a SchurHelper from a given DomainCollection
	 *
	 * @param domain the DomainCollection
	 * @param comm the teuchos communicator
	 */
	SchurHelper(std::shared_ptr<Domain<D>> domain, std::shared_ptr<PatchSolver<D>> solver,
	            std::shared_ptr<PatchOperator<D>> op, std::shared_ptr<IfaceInterp<D>> interpolator);

	/**
	 * @brief Solve with a given set of interface values
	 *
	 * @param f the rhs vector
	 * @param u the vector to put solution in
	 * @param gamma the interface values to use
	 * @param diff the resulting difference
	 */
	void solveWithInterface(std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u,
	                        std::shared_ptr<const Vector<D - 1>> gamma,
	                        std::shared_ptr<Vector<D - 1>>       diff);
	void solveAndInterpolateWithInterface(std::shared_ptr<const Vector<D>>     f,
	                                      std::shared_ptr<Vector<D>>           u,
	                                      std::shared_ptr<const Vector<D - 1>> gamma,
	                                      std::shared_ptr<Vector<D - 1>>       interp);
	void solveWithSolution(std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u);
	void interpolateToInterface(std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u,
	                            std::shared_ptr<Vector<D - 1>> gamma);

	/**
	 * @brief Apply patch operator with a given set of interface values
	 *
	 * @param u the solution vector to use
	 * @param gamma the interface values to use
	 * @param f the resulting rhs vector
	 */
	void applyWithInterface(std::shared_ptr<const Vector<D>>     u,
	                        std::shared_ptr<const Vector<D - 1>> gamma,
	                        std::shared_ptr<Vector<D>>           f);
	void apply(std::shared_ptr<const Vector<D>> u, std::shared_ptr<Vector<D>> f);

	void scatterInterface(std::shared_ptr<Vector<D - 1>>       gamma_dist,
	                      std::shared_ptr<const Vector<D - 1>> gamma)
	{
		// TODO make this general;
		PetscVector<D - 1> *      p_dist   = dynamic_cast<PetscVector<D - 1> *>(gamma_dist.get());
		const PetscVector<D - 1> *p_global = dynamic_cast<const PetscVector<D - 1> *>(gamma.get());
		if (p_dist == nullptr || p_global == nullptr) { throw 3; }
		VecScatterBegin(scatter, p_global->vec, p_dist->vec, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(scatter, p_global->vec, p_dist->vec, INSERT_VALUES, SCATTER_FORWARD);
	}

	void scatterInterfaceReverse(std::shared_ptr<const Vector<D - 1>> gamma_dist,
	                             std::shared_ptr<Vector<D - 1>>       gamma)
	{
		// TODO make this general;
		const PetscVector<D - 1> *p_dist
		= dynamic_cast<const PetscVector<D - 1> *>(gamma_dist.get());
		PetscVector<D - 1> *p_global = dynamic_cast<PetscVector<D - 1> *>(gamma.get());
		if (p_dist == nullptr || p_global == nullptr) { throw 3; }
		VecScatterBegin(scatter, p_dist->vec, p_global->vec, ADD_VALUES, SCATTER_REVERSE);
		VecScatterEnd(scatter, p_dist->vec, p_global->vec, ADD_VALUES, SCATTER_REVERSE);
	}
	void updateInterfaceDist(std::shared_ptr<Vector<D - 1>> gamma_dist)
	{
		gamma->set(0);
		scatterInterfaceReverse(gamma_dist, gamma);
		scatterInterface(gamma_dist, gamma);
	}
	std::shared_ptr<PetscVector<D - 1>> getNewSchurVec()
	{
		Vec u;
		VecCreateMPI(MPI_COMM_WORLD, iface_map_vec.size() * iface_stride, PETSC_DETERMINE, &u);
		return std::shared_ptr<PetscVector<D - 1>>(new PetscVector<D - 1>(u, lengths));
	}
	std::shared_ptr<PetscVector<D - 1>> getNewSchurDistVec()
	{
		Vec u;
		VecCreateSeq(PETSC_COMM_SELF, iface_dist_map_vec.size() * iface_stride, &u);
		return std::shared_ptr<PetscVector<D - 1>>(new PetscVector<D - 1>(u, lengths));
	}

	int getSchurVecLocalSize() const
	{
		return iface_map_vec.size() * iface_stride;
	}
	int getSchurVecGlobalSize() const
	{
		return num_global_ifaces * iface_stride;
	}
	// getters
	std::shared_ptr<IfaceInterp<D>> getIfaceInterp()
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
	const std::array<int, D - 1> getLengths() const
	{
		return lengths;
	}
};
template <size_t D>
inline SchurHelper<D>::SchurHelper(std::shared_ptr<Domain<D>>        domain,
                                   std::shared_ptr<PatchSolver<D>>   solver,
                                   std::shared_ptr<PatchOperator<D>> op,
                                   std::shared_ptr<IfaceInterp<D>>   interpolator)
{
	iface_stride = 1;
	for (size_t i = 0; i < D - 1; i++) {
		iface_stride *= domain->getNs()[i];
		lengths[i] = domain->getNs()[i];
	}
	sinfo_vector.reserve(domain->getNumLocalPatches());
	for (auto &pinfo : domain->getPatchInfoVector()) {
		sinfo_vector.push_back(pinfo);
	}
	std::map<int, std::map<int, IfaceSet<D>>> off_proc_ifaces;
	std::set<int>                             incoming_procs;
	for (SchurInfo<D> &sd : sinfo_vector) {
		sd.enumerateIfaces(ifaces, off_proc_ifaces, incoming_procs);
		solver->addDomain(sd);
	}
	{
		using namespace std;
		// send info
		deque<char *>       buffers;
		vector<MPI_Request> send_requests;
		for (auto &p : off_proc_ifaces) {
			int dest = p.first;
			int size = 0;
			for (auto q : p.second) {
				IfaceSet<D> &iface = q.second;
				size += iface.serialize(nullptr);
			}
			char *buffer = new char[size];
			buffers.push_back(buffer);
			int pos = 0;
			for (auto q : p.second) {
				IfaceSet<D> &iface = q.second;
				pos += iface.serialize(buffer + pos);
			}
			MPI_Request request;
			MPI_Isend(buffer, size, MPI_BYTE, dest, 0, MPI_COMM_WORLD, &request);
			send_requests.push_back(request);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// recv info
		for (int src : incoming_procs) {
			MPI_Status status;
			MPI_Probe(src, 0, MPI_COMM_WORLD, &status);
			int size;
			MPI_Get_count(&status, MPI_BYTE, &size);
			char *buffer = new char[size];

			MPI_Recv(buffer, size, MPI_BYTE, src, 0, MPI_COMM_WORLD, &status);

			BufferReader reader(buffer);
			while (reader.getPos() < size) {
				IfaceSet<D> ifs;
				reader >> ifs;
				ifaces[ifs.id].insert(ifs);
			}

			delete[] buffer;
		}
		// wait for all
		MPI_Waitall(send_requests.size(), &send_requests[0], MPI_STATUSES_IGNORE);
		// delete send buffers
		for (char *buffer : buffers) {
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
	gamma              = getNewSchurVec();
	PW<IS> dist_is;
	ISCreateBlock(MPI_COMM_SELF, iface_stride, iface_dist_map_vec.size(), &iface_dist_map_vec[0],
	              PETSC_COPY_VALUES, &dist_is);
	VecScatterCreate(gamma->vec, dist_is, local_gamma->vec, nullptr, &scatter);

	int num_ifaces = ifaces.size();
	MPI_Allreduce(&num_ifaces, &num_global_ifaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}
template <size_t D>
inline void SchurHelper<D>::solveWithInterface(std::shared_ptr<const Vector<D>>     f,
                                               std::shared_ptr<Vector<D>>           u,
                                               std::shared_ptr<const Vector<D - 1>> gamma,
                                               std::shared_ptr<Vector<D - 1>>       diff)
{
	scatterInterface(local_gamma, gamma);

	solver->domainSolve(sinfo_vector, f, u, local_gamma);

	local_gamma->set(0);
	for (SchurInfo<D> &sd : sinfo_vector) {
		interpolator->interpolate(sd, u, local_gamma);
	}

	diff->set(0);
	scatterInterfaceReverse(local_gamma, diff);
	diff->addScaled(-1, gamma);
}
template <size_t D>
inline void SchurHelper<D>::solveAndInterpolateWithInterface(
std::shared_ptr<const Vector<D>> f, std::shared_ptr<Vector<D>> u,
std::shared_ptr<const Vector<D - 1>> gamma, std::shared_ptr<Vector<D - 1>> interp)
{
	scatterInterface(local_gamma, gamma);

	// solve over sinfo_vector on this proc
	solver->domainSolve(sinfo_vector, f, u, local_gamma);

	local_gamma->set(0);
	for (SchurInfo<D> &sd : sinfo_vector) {
		interpolator->interpolate(sd, u, local_gamma);
	}

	interp->set(0);
	scatterInterfaceReverse(local_gamma, interp);
}
template <size_t D>
inline void SchurHelper<D>::solveWithSolution(std::shared_ptr<const Vector<D>> f,
                                              std::shared_ptr<Vector<D>>       u)
{
	local_gamma->set(0);
	for (SchurInfo<D> &sd : sinfo_vector) {
		interpolator->interpolate(sd, u, local_gamma);
	}

	updateInterfaceDist(local_gamma);

	// solve over sinfo_vector on this proc
	solver->domainSolve(sinfo_vector, f, u, local_gamma);
}
template <size_t D>
inline void SchurHelper<D>::interpolateToInterface(std::shared_ptr<const Vector<D>> f,
                                                   std::shared_ptr<Vector<D>>       u,
                                                   std::shared_ptr<Vector<D - 1>>   gamma)
{
	// initilize our local variables
	local_gamma->set(0);
	for (SchurInfo<D> &sd : sinfo_vector) {
		interpolator->interpolate(sd, u, local_gamma);
	}
	gamma->set(0);
	scatterInterfaceReverse(local_gamma, gamma);
}
template <size_t D>
inline void SchurHelper<D>::applyWithInterface(std::shared_ptr<const Vector<D>>     u,
                                               std::shared_ptr<const Vector<D - 1>> gamma,
                                               std::shared_ptr<Vector<D>>           f)
{
	scatterInterface(local_gamma, gamma);

	f->set(0);

	for (SchurInfo<D> &sd : sinfo_vector) {
		op->apply(sd, u, local_gamma, f);
	}
}
template <size_t D>
inline void SchurHelper<D>::apply(std::shared_ptr<const Vector<D>> u, std::shared_ptr<Vector<D>> f)
{
	local_gamma->set(0);
	for (SchurInfo<D> &sd : sinfo_vector) {
		interpolator->interpolate(sd, u, local_gamma);
	}

	updateInterfaceDist(local_gamma);

	f->set(0);
	for (SchurInfo<D> &sd : sinfo_vector) {
		op->apply(sd, u, local_gamma, f);
	}
}
template <size_t D> inline void SchurHelper<D>::indexDomainIfacesLocal()
{
	using namespace std;
	vector<int>   map_vec;
	map<int, int> rev_map;
	if (!sinfo_vector.empty()) {
		int curr_i = 0;
		for (SchurInfo<D> &sd : sinfo_vector) {
			for (int id : sd.getIds()) {
				if (rev_map.count(id) == 0) {
					rev_map[id] = curr_i;
					map_vec.push_back(id);
					curr_i++;
				}
			}
		}
		for (SchurInfo<D> &sd : sinfo_vector) {
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
		for (SchurInfo<D> &sd : sinfo_vector) {
			sd.setGlobalIndexes(rev_map);
		}
		for (size_t i = 0; i < iface_dist_map_vec.size(); i++) {
			iface_dist_map_vec[i] = inds[i];
		}
	}
}
template <size_t D> class SchurHelperVG : public VectorGenerator<D>
{
	private:
	std::shared_ptr<SchurHelper<D + 1>> sh;

	public:
	SchurHelperVG(std::shared_ptr<SchurHelper<D + 1>> sh)
	{
		this->sh = sh;
	}
	std::shared_ptr<Vector<D>> getNewVector()
	{
		return sh->getNewSchurVec();
	}
};
extern template class SchurHelper<2>;
extern template class SchurHelper<3>;
#endif
