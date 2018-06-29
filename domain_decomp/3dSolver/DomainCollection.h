#ifndef DOMAINSIGNATURECOLLECTION_H
#define DOMAINSIGNATURECOLLECTION_H
#include "Domain.h"
#include "InterpCase.h"
#include "OctTree.h"
#include "PW.h"
#include "Side.h"
#include <map>
#include <memory>
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
class DomainCollection
{
	private:
	int  n = -1;
	void indexDomainsLocal();
	void indexDomainsGlobal();
	void zoltanBalanceDomains();
	void reIndex();

	public:
	int  rank;
	int  num_pins;
	bool neumann = false;
	/**
	 * @brief Number of total domains.
	 */
	int num_global_domains = 1;
	/**
	 * @brief Number of total interfaces.
	 */
	int num_global_interfaces = 0;
	/**
	 * @brief A map that maps the id of a domain to its domain signature.
	 */
	std::map<int, std::shared_ptr<Domain>> domains;

	std::vector<int> domain_gid_map_vec;
	std::vector<int> domain_map_vec;
	std::vector<int> domain_off_proc_map_vec;

	int getN() { return n; }
	/**
	 * @brief Default empty constructor.
	 */
	DomainCollection() = default;

	DomainCollection(OctTree t, int n);
	DomainCollection(OctTree t, int level, int n);
	/**
	 * @brief Generate a grid of domains.
	 *
	 * @param d_x number of domains in the x direction.
	 * @param d_y number of domains in the y direction.
	 * @param rank the rank of the MPI process.
	 */
	DomainCollection(int d_x, int d_y, int d_z, int n);
	/**
	 * @brief Balance the domains over processors using Zoltan
	 */
	void zoltanBalance();
	void zoltanBalanceWithLower(DomainCollection &lower);


	void setNeumann()
	{
		neumann = true;
		for (auto &p : domains) {
			p.second->setNeumann();
		}
	}
	PW_explicit<Vec> getNewDomainVec() const;

	int    getGlobalNumCells() { return num_global_domains * n * n * n; }
	double integrate(const Vec u);
	double volume();
};
#endif
