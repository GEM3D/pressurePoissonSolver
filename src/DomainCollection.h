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

#ifndef DOMAINSIGNATURECOLLECTION_H
#define DOMAINSIGNATURECOLLECTION_H
#include <Domain.h>
#include <InterpCase.h>
#include <OctTree.h>
#include <PW.h>
#include <PetscVector.h>
#include <cmath>
#include <deque>
#include <map>
#include <memory>
#include <petscao.h>
#include <petscvec.h>
#include <set>
#include <string>
#include <vector>
/**
 * @brief A collection of Domain Signatures
 *
 * There are two purposes for this class:
 *   -# Partition the domains across processors
 *   -# Determine the number of interfaces, and provide a unique index to each interface.
 */
template <size_t D> class DomainCollection
{
	private:
	std::array<int, D> lengths;
	int                patch_stride;
	void               indexDomainsLocal()
	{
		int                curr_i = 0;
		std::vector<int>   map_vec;
		std::vector<int>   off_proc_map_vec;
		std::map<int, int> rev_map;
		std::set<int>      offs;
		if (!domains.empty()) {
			std::set<int> todo;
			for (auto &p : domains) {
				todo.insert(p.first);
			}
			std::set<int> enqueued;
			while (!todo.empty()) {
				std::deque<int> queue;
				queue.push_back(*todo.begin());
				enqueued.insert(*todo.begin());
				while (!queue.empty()) {
					int i = queue.front();
					todo.erase(i);
					queue.pop_front();
					map_vec.push_back(i);
					Domain<D> &d = *domains[i];
					rev_map[i]   = curr_i;
					d.id_local   = curr_i;
					curr_i++;
					for (int i : d.getNbrIds()) {
						if (!enqueued.count(i)) {
							enqueued.insert(i);
							if (domains.count(i)) {
								queue.push_back(i);
							} else {
								if (!offs.count(i)) {
									offs.insert(i);
									off_proc_map_vec.push_back(i);
								}
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
		for (auto &p : domains) {
			p.second->setLocalNeighborIndexes(rev_map);
		}
		// domain_rev_map          = rev_map;
		domain_map_vec          = map_vec;
		domain_gid_map_vec      = map_vec;
		domain_off_proc_map_vec = off_proc_map_vec;
		indexDomainsGlobal();
	}
	void indexDomainsGlobal()
	{
		// global indices are going to be sequentially increasing with rank
		int local_size = domains.size();
		int start_i;
		MPI_Scan(&local_size, &start_i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		start_i -= local_size;
		std::vector<int> new_global(local_size);
		iota(new_global.begin(), new_global.end(), start_i);

		// create map for gids
		PW<AO> ao;
		AOCreateMapping(MPI_COMM_WORLD, local_size, &domain_map_vec[0], &new_global[0], &ao);

		// get global indices that we want to recieve for dest vector
		std::vector<int> inds = domain_map_vec;
		for (int i : domain_off_proc_map_vec) {
			inds.push_back(i);
		}

		// get new global indices
		AOApplicationToPetsc(ao, inds.size(), &inds[0]);
		std::map<int, int> rev_map;
		for (size_t i = 0; i < inds.size(); i++) {
			rev_map[i] = inds[i];
		}

		for (auto &p : domains) {
			p.second->setGlobalNeighborIndexes(rev_map);
		}
		for (size_t i = 0; i < domain_map_vec.size(); i++) {
			domain_map_vec[i] = inds[i];
		}
		for (size_t i = 0; i < domain_off_proc_map_vec.size(); i++) {
			domain_off_proc_map_vec[i] = inds[domain_map_vec.size() + i];
		}
	}

	void zoltanBalanceDomains();
	void reIndex()
	{
		indexDomainsLocal();
	}

	public:
	int  rank;
	bool neumann = false;
	/**
	 * @brief Number of total domains.
	 */
	int num_global_domains = 1;
	/**
	 * @brief A map that maps the id of a domain to its domain signature.
	 */
	std::map<int, std::shared_ptr<Domain<D>>>        domains;
	const std::map<int, std::shared_ptr<Domain<D>>> &getDomainMap()
	{
		return domains;
	}

	std::vector<int> domain_gid_map_vec;
	std::vector<int> domain_map_vec;
	std::vector<int> domain_off_proc_map_vec;

	const std::array<int, D> &getLengths()
	{
		return lengths;
	}
	/**
	 * @brief Default empty constructor.
	 */
	DomainCollection() = default;
	DomainCollection(std::map<int, std::shared_ptr<Domain<D>>> domain_set,
	                 const std::array<int, D> &                lengths)
	{
		this->lengths = lengths;

		patch_stride = 1;
		for (int i = 0; i < D; i++) {
			patch_stride *= lengths[i];
		}

		domains = domain_set;

		int num_local_domains = domains.size();
		MPI_Allreduce(&num_local_domains, &num_global_domains, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		reIndex();

		for (auto &p : domains) {
			p.second->setPtrs(domains);
		}
	}

	void setNeumann()
	{
		neumann = true;
		for (auto &p : domains) {
			p.second->setNeumann();
		}
	}
	std::shared_ptr<PetscVector<D>> getNewDomainVec() const
	{
		Vec u;
		VecCreateMPI(MPI_COMM_WORLD, domains.size() * patch_stride, PETSC_DETERMINE, &u);
		return std::shared_ptr<PetscVector<D>>(new PetscVector<D>(u, lengths));
	}

	int getGlobalNumDomains()
	{
		return num_global_domains;
	}
	int getGlobalNumCells()
	{
		return num_global_domains * patch_stride;
	}
	int getLocalNumCells()
	{
		return domains.size() * patch_stride;
	}
	int getNumElementsInDomain()
	{
		return patch_stride;
	}
	double volume()
	{
		double sum = 0;
		for (auto &p : domains) {
			Domain<D> &d = *p.second;
			sum
			+= std::accumulate(d.lengths.begin(), d.lengths.end(), 1.0, std::multiplies<double>());
		}
		double retval;
		MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return retval;
	}
	double integrate(std::shared_ptr<const Vector<D>> u) const
	{
		double sum = 0;

		for (auto &p : domains) {
			Domain<D> &        d      = *p.second;
			const LocalData<D> u_data = u->getLocalData(d.id_local);

			double patch_sum = 0;
			nested_loop<D>(u_data.getStart(), u_data.getEnd(),
			               [&](std::array<int, D> coord) { patch_sum += u_data[coord]; });

			patch_sum
			*= std::accumulate(d.lengths.begin(), d.lengths.end(), 1.0, std::multiplies<double>())
			   / std::pow(d.n, D);

			sum += patch_sum;
		}
		double retval;
		MPI_Allreduce(&sum, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		return retval;
	}
};
#endif
