#include "DomainCollection.h"
#include <tuple>
#include <array>
using Teuchos::RCP;
using Teuchos::rcp;
using namespace std;
extern "C" {
// LU decomoposition of a general matrix
void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);

// generate inverse of a matrix given its LU decomposition
void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork, int *INFO);
}

enum axis_enum { X_AXIS, Y_AXIS };
enum bc_enum { DIRICHLET, NEUMANN, REFINED };

class Iface
{
	public:
	bool        right;
	int         axis;
    array<int,8> global_i;
    array<int,4> types;
	friend bool operator<(const Iface &l, const Iface &r)
	{
		return std::tie(l.global_i[0], l.right) < std::tie(r.global_i[0], r.right);
	}
	friend bool operator==(const Iface &l, const Iface &r)
	{
		return l.types == r.types;
	}
	/*friend bool operator!=(const Iface &l, const Iface &r)
	{
		return std::tie(l.l_south, l.t_south, l.l_east, l.t_east, l.l_north, l.t_north, l.l_west,
		                l.t_west)
		       == std::tie(r.l_south, r.t_south, r.l_west, r.t_west, r.l_north, r.t_north, r.l_east,
		                   r.t_east);
	}*/
};

DomainCollection::DomainCollection(DomainSignatureCollection dsc, int n, double h_x, double h_y,
                                   RCP<const Teuchos::Comm<int>> comm)
{
	// cerr<< "Low:  " << low << "\n";
	// cerr<< "High: " << high << "\n";
	this->comm         = comm;
	this->n           = n;
	this->h_x          = h_x;
	this->h_y          = h_y;
	num_global_domains = dsc.num_global_domains;
	for (auto p : dsc.domains) {
		DomainSignature ds       = p.second;
		int             i        = ds.id;

		// create a domain
		RCP<Domain> d_ptr = rcp(new Domain(ds, n, n, h_x, h_y));
		domains[i]        = d_ptr;
	}
}

void DomainCollection::initNeumann(function<double(double, double)> ffun,
                                   function<double(double, double)> efun,
                                   function<double(double, double)> nfunx,
                                   function<double(double, double)> nfuny)
{
	for (auto &p : domains) {
		Domain &d        = *p.second;

		// Generate RHS vector
		std::valarray<double> &f     = d.f;
		std::valarray<double> &exact = d.exact;

		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x           = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = efun(x, y);
				// west
				if (d.nbr[6] == -1) {
					d.boundary_west[yi] = nfunx(0.0, y);
				}
				// east
				if (d.nbr[2] == -1) {
					d.boundary_east[yi] = nfunx(1.0, y);
				}
				// south
				if (d.nbr[4] == -1) {
					d.boundary_south[xi] = nfuny(x, 0.0);
				}
				// north
				if (d.nbr[0] == -1) {
					d.boundary_north[xi] = nfuny(x, 1.0);
				}
			}
		}
		d.planNeumann();
	}

	// create map for domains
	generateMaps();
	distributeIfaceInfo();
}

void DomainCollection::initDirichlet(function<double(double, double)> ffun,
                                     function<double(double, double)> gfun)
{
	for (auto &p : domains) {
		Domain &d        = *p.second;

		// Generate RHS vector
		std::valarray<double> &f     = d.f;
		std::valarray<double> &exact = d.exact;
		// Use local indices to access the entries of f_data.
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x           = d.x_start + h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = gfun(x, y);
				if (d.nbr[6] == -1) {
					d.boundary_west[yi] = gfun(0.0, y);
				}
				if (d.nbr[2] == -1) {
					d.boundary_east[yi] = gfun(1.0, y);
				}
				if (d.nbr[4] == -1) {
					d.boundary_south[xi] = gfun(x, 0.0);
				}
				if (d.nbr[0] == -1) {
					d.boundary_north[xi] = gfun(x, 1.0);
				}
			}
		}
		d.planDirichlet();
	}
	// create map for domains
	generateMaps();
	distributeIfaceInfo();
}
void DomainCollection::initDirichletRefined(function<double(double, double)> ffun,
                                     function<double(double, double)> gfun)
{
	amr = true;
	for (auto &p : domains) {
		Domain &d = *p.second;

		// Generate RHS vector
		std::valarray<double> &f     = d.f;
		std::valarray<double> &exact = d.exact;
		// Use local indices to access the entries of f_data.
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x           = d.x_start + d.h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + d.h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = gfun(x, y);
				if (d.nbr[6] == -1) {
					d.boundary_west[yi] = gfun(0.0, y);
				}
				if (d.nbr[2] == -1) {
					d.boundary_east[yi] = gfun(2.0, y);
				}
				if (d.nbr[4] == -1) {
					d.boundary_south[xi] = gfun(x, 0.0);
				}
				if (d.nbr[0] == -1) {
					d.boundary_north[xi] = gfun(x, 1.0);
				}
			}
		}
		d.planDirichlet();
	}
	// create map for domains
	generateMaps();
	//distributeIfaceInfo();
}
void DomainCollection::generateMaps()
{
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	set<int>   not_visited;
	for (auto &p : domains) {
		not_visited.insert(p.first);
	}
	vector<int> global;
	int         curr_i = 0;
	vector<int> c_iface_global;
	int         curr_c_i = 0;
	vector<int> matrix_global;
	int         curr_matrix_i = 0;
	vector<int> iface_global;
	global.reserve(domains.size() * (2 * n + 2 * n));
	while (!not_visited.empty()) {
		int first = *not_visited.begin();
		queue.push_back(first);
		enqueued.insert(first);
		while (!queue.empty()) {
			int     curr = queue.front();
			Domain &d    = *domains[curr];
			queue.pop_front();
			visited.insert(curr);
            not_visited.erase(curr);
			for (int q = 0; q < 8; q++) {
				if (d.nbr[q] != -1 && visited.count(d.nbr[q]) == 0) {
					// a new edge that we have not assigned an index to
					try {
						Domain &nbr   = *domains.at(d.nbr[q]);
						int     nbr_q = -1;
						do {
							nbr_q++;
						} while (nbr.nbr[nbr_q] != curr);
						nbr.local_i[nbr_q]       = curr_i;
						nbr.iface_local_i[nbr_q] = curr_c_i;
						if (enqueued.count(d.nbr[q]) == 0) {
							queue.push_back(d.nbr[q]);
							enqueued.insert(d.nbr[q]);
						}
					} catch (const out_of_range &oor) {
						// do nothing
					}
					d.local_i[q]       = curr_i;
					d.iface_local_i[q] = curr_c_i;
					for (int i = 0; i < 22; i++) {
						c_iface_global.push_back(d.iface_i[q] + i);
						iface_global.push_back(d.iface_i[q] + i);
						curr_c_i++;
					}
					for (int i = 0; i < n; i++) {
						global.push_back(d.global_i[q] + i);
						matrix_global.push_back(d.global_i[q] + i);
						curr_matrix_i++;
					}
					curr_i += n;
				}
			}
		}
	}
	// Now that the global indices have been calculated, we can create a map for the interface
	// points
	if (num_global_domains == 1) {
		// this is a special case for when there is only one domain
		collection_map       = Teuchos::rcp(new map_type(1, 0, comm));
		collection_iface_map = Teuchos::rcp(new map_type(1, 0, comm));
		matrix_map           = Teuchos::rcp(new map_type(1, 0, comm));
		iface_map            = Teuchos::rcp(new map_type(1, 0, comm));
	} else {
		collection_map       = Teuchos::rcp(new map_type(-1, &global[0], global.size(), 0, this->comm));
		collection_iface_map
		= Teuchos::rcp(new map_type(-1, &c_iface_global[0], c_iface_global.size(), 0, this->comm));
		matrix_map
		= Teuchos::rcp(new map_type(-1, &matrix_global[0], matrix_global.size(), 0, this->comm));
		iface_map
		= Teuchos::rcp(new map_type(-1, &iface_global[0], iface_global.size(), 0, this->comm));
	}
}
void DomainCollection::distributeIfaceInfo(){
    //
	int_vector_type dist(collection_iface_map, 1);
	iface_info = rcp(new int_vector_type(iface_map, 1));
	auto             dist_view = dist.getLocalView<Kokkos::HostSpace>();
	for (auto &p : domains) {
		Domain &d2 = *p.second;
		if (d2.nbr[0] != -1) {
			//dist_view(d2.iface_local_i[0], 0)      = d2.global_i[0];
			//dist_view(d2.iface_local_i[0] + 1, 0)  = 0;
			//dist_view(d2.iface_local_i[0] + 2, 0)  = n;
			dist_view(d2.iface_local_i[0] + 12, 0) = d2.global_i[2];
			if (d2.neumann && d2.nbr[2] == -1) {
				dist_view(d2.iface_local_i[0] + 13, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[0] + 13, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[0] + 14, 0) = n;
			dist_view(d2.iface_local_i[0] + 15, 0) = d2.global_i[4];
			if (d2.neumann && d2.nbr[4] == -1) {
				dist_view(d2.iface_local_i[0] + 16, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[0] + 16, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[0] + 17, 0) = n;
			dist_view(d2.iface_local_i[0] + 18, 0) = d2.global_i[6];
			if (d2.neumann && d2.nbr[6] == -1) {
				dist_view(d2.iface_local_i[0] + 19, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[0] + 19, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[0] + 20, 0) = n;
			dist_view(d2.iface_local_i[0] + 21, 0) = X_AXIS;
		}
		if (d2.nbr[2] != -1) {
			//dist_view(d2.iface_local_i[2], 0)      = d2.global_i[2];
			//dist_view(d2.iface_local_i[2] + 1, 0)  = 0;
			//dist_view(d2.iface_local_i[2] + 2, 0)  = n;
			dist_view(d2.iface_local_i[2] + 12, 0) = d2.global_i[4];
			if (d2.neumann && d2.nbr[4] == -1) {
				dist_view(d2.iface_local_i[2] + 13, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[2] + 13, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[2] + 14, 0) = n;
			dist_view(d2.iface_local_i[2] + 15, 0) = d2.global_i[6];
			if (d2.neumann && d2.nbr[6] == -1) {
				dist_view(d2.iface_local_i[2] + 16, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[2] + 16, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[2] + 17, 0) = n;
			dist_view(d2.iface_local_i[2] + 18, 0) = d2.global_i[0];
			if (d2.neumann && d2.nbr[0] == -1) {
				dist_view(d2.iface_local_i[2] + 19, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[2] + 19, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[2] + 20, 0) = n;
			dist_view(d2.iface_local_i[2] + 21, 0) = Y_AXIS;
		}
		if (d2.nbr[4] != -1) {
			dist_view(d2.iface_local_i[4], 0)      = d2.global_i[4];
			dist_view(d2.iface_local_i[4] + 1, 0)  = 0;
			if (d2.neumann && d2.nbr[4] == -1) {
				dist_view(d2.iface_local_i[4] + 1, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[4] + 1, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[4] + 2, 0)  = n;
			dist_view(d2.iface_local_i[4] + 3, 0)  = d2.global_i[6];
			if (d2.neumann && d2.nbr[6] == -1) {
				dist_view(d2.iface_local_i[4] + 4, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[4] + 4, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[4] + 5, 0)  = n;
			dist_view(d2.iface_local_i[4] + 6, 0)  = d2.global_i[0];
			if (d2.neumann && d2.nbr[0] == -1) {
				dist_view(d2.iface_local_i[4] + 7, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[4] + 7, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[4] + 8, 0)  = n;
			dist_view(d2.iface_local_i[4] + 9, 0)  = d2.global_i[2];
			if (d2.neumann && d2.nbr[2] == -1) {
				dist_view(d2.iface_local_i[4] + 10, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[4] + 10, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[4] + 11, 0) = n;
			//dist_view(d2.iface_local_i[4] + 21, 0) = X_AXIS;
		}
		if (d2.nbr[6] != -1) {
			dist_view(d2.iface_local_i[6], 0)      = d2.global_i[6];
			if (d2.neumann && d2.nbr[6] == -1) {
				dist_view(d2.iface_local_i[6] + 1, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[6] + 1, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[6] + 2, 0) = n;
			dist_view(d2.iface_local_i[6] + 3, 0) = d2.global_i[0];
			if (d2.neumann && d2.nbr[0] == -1) {
				dist_view(d2.iface_local_i[6] + 4, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[6] + 4, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[6] + 5, 0) = n;
			dist_view(d2.iface_local_i[6] + 6, 0) = d2.global_i[2];
			if (d2.neumann && d2.nbr[2] == -1) {
				dist_view(d2.iface_local_i[6] + 7, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[6] + 7, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[6] + 8, 0) = n;
			dist_view(d2.iface_local_i[6] + 9, 0) = d2.global_i[4];
			if (d2.neumann && d2.nbr[4] == -1) {
				dist_view(d2.iface_local_i[6] + 10, 0) = NEUMANN;
			} else {
				dist_view(d2.iface_local_i[6] + 10, 0) = DIRICHLET;
			}
			dist_view(d2.iface_local_i[6] + 11, 0) = n;
			// dist_view(d2.iface_local_i[6] + 21, 0) = Y_AXIS;
		}
	}
	Tpetra::Export<> exporter(collection_iface_map, iface_map);
	iface_info->doExport(dist, exporter, Tpetra::CombineMode::ADD);
}
void DomainCollection::solveWithInterface(const vector_type &gamma, vector_type &diff)
{
	// auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
	// if(has_east)std::cout << "Gamma begin\n";
	// gamma.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	diff.update(1, gamma, 0);
	Tpetra::Import<> importer(diff.getMap(), collection_map);
	Tpetra::Export<> exporter(collection_map, diff.getMap());
	vector_type      local_gamma(collection_map, 1);
	vector_type      local_diff(collection_map, 1);
	local_gamma.doImport(diff, importer, Tpetra::CombineMode::INSERT);

	// solve over domains on this proc
	for (auto &p : domains) {
		p.second->solveWithInterface(local_gamma, local_diff);
	}

	// export diff vector

	// if(has_east)std::cout <<"LOCAL AFEr\n";
	// local_vector.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	diff.scale(0);
	// if(has_east)std::cout <<"DIFFBEFORE\n";
	// diff.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	diff.doExport(local_diff, importer, Tpetra::CombineMode::ADD);
	// if(has_east)std::cout <<"DIFF\n";
	// diff.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	// if(has_east)std::cout <<"GAMMA\n";
	// gamma.describe(*out,Teuchos::EVerbosityLevel::VERB_EXTREME);
	//diff.update(-2, gamma, 1);
}
double DomainCollection::diffNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second->diffNorm(), 2);
	}
	return sqrt(result);
}
double DomainCollection::diffNorm(double uavg, double eavg)
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second->diffNorm(uavg, eavg), 2);
	}
	return sqrt(result);
}
double DomainCollection::uSum()
{
	double result = 0;
	for (auto &p : domains) {
		result += p.second->uSum();
	}
	return result;
}
double DomainCollection::exactNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second->exactNorm(), 2);
	}
	return sqrt(result);
}
double DomainCollection::fNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second->fNorm(), 2);
	}
    double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &result, &retval);
	return sqrt(retval);
}
double DomainCollection::exactNorm(double eavg)
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second->exactNorm(eavg), 2);
	}
	return sqrt(result);
}
double DomainCollection::exactSum()
{
	double result = 0;
	for (auto &p : domains) {
		result += p.second->exactSum();
	}
	return result;
}
double DomainCollection::residual()
{
	vector_type      ghost(collection_map, 2);
	vector_type      one_ghost(matrix_map, 2);
	for (auto &p : domains) {
		p.second->putGhostCells(ghost);
	}
	Tpetra::Export<> exporter(collection_map, matrix_map);
	one_ghost.doExport(ghost, exporter, Tpetra::CombineMode::ADD);
	ghost.putScalar(0);
	ghost.doImport(one_ghost, exporter, Tpetra::CombineMode::ADD);

	double residual = 0;
	for (auto &p : domains) {
		residual += pow(p.second->residual(ghost), 2);
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &residual, &retval);
	return sqrt(retval);
}
RCP<matrix_type> DomainCollection::formMatrix(RCP<map_type> map, int delete_row)
{
	int              size = max(n, n);
	RCP<matrix_type> A    = rcp(new matrix_type(map, size * 6));
    // create iface objects
	/*set<Iface> ifaces;
	auto       iface_view = iface_info->getLocalView<Kokkos::HostSpace>();
	for (size_t i = 0; i < iface_view.dimension(0); i += 22) {
		Iface right;
		Iface left;

		right.right   = true;
		right.axis    = iface_view(i + 21, 0);
		right.i_south = iface_view(i, 0);
		right.t_south = iface_view(i + 1, 0);
		right.l_south = iface_view(i + 2, 0);
		right.i_west  = iface_view(i + 3, 0);
		right.t_west  = iface_view(i + 4, 0);
		right.l_west  = iface_view(i + 5, 0);
		right.i_north = iface_view(i + 6, 0);
		right.t_north = iface_view(i + 7, 0);
		right.l_north = iface_view(i + 8, 0);
		right.i_east  = iface_view(i + 9, 0);
		right.t_east  = iface_view(i + 10, 0);
		right.l_east  = iface_view(i + 11, 0);

		left.right   = false;
		left.axis    = iface_view(i + 21, 0);
		left.i_south = iface_view(i, 0);
		left.t_south = iface_view(i + 1, 0);
		left.l_south = iface_view(i + 2, 0);
		left.i_west  = iface_view(i + 12, 0);
		left.t_west  = iface_view(i + 13, 0);
		left.l_west  = iface_view(i + 14, 0);
		left.i_north = iface_view(i + 15, 0);
		left.t_north = iface_view(i + 16, 0);
		left.l_north = iface_view(i + 17, 0);
		left.i_east  = iface_view(i + 18, 0);
		left.t_east  = iface_view(i + 19, 0);
		left.l_east  = iface_view(i + 20, 0);
        
        ifaces.insert(left);
        ifaces.insert(right);
	}
    
    int num_types=0;
	while(!ifaces.empty()){
        num_types++;
		// the first in the set is the type of interface that we are going to solve for
		set<Iface> todo;
		Iface      curr_type = *ifaces.begin();
		ifaces.erase(ifaces.begin());
		todo.insert(curr_type);
        set<Iface> to_be_deleted;
		for (auto iter = ifaces.begin(); iter != ifaces.end(); iter++) {
			if (*iter == curr_type) {
				todo.insert(*iter);
				to_be_deleted.insert(*iter);

				// TODO fix this iterator
				// iter=ifaces.begin();
			}
		}
		for (Iface i : to_be_deleted) {
			ifaces.erase(i);
		}

		// create domain representing curr_type
		DomainSignature ds;
		if (curr_type.t_north == NEUMANN) {
			ds.nbr[0] = -1;
		} else {
			ds.nbr[0] = 1;
		}
		if (curr_type.t_east == NEUMANN) {
			ds.nbr[2] = -1;
		} else {
			ds.nbr[2] = 1;
		}
		if (curr_type.t_south == NEUMANN) {
			ds.nbr[4] = -1;
		} else {
			ds.nbr[4] = 1;
		}
		if (curr_type.t_west == NEUMANN) {
			ds.nbr[6] = -1;
		} else {
			ds.nbr[6] = 1;
		}
		Domain d(ds, n, n, h_x, h_y);

		d.boundary_north = valarray<double>(n);
		d.boundary_south = valarray<double>(n);
		d.boundary_east  = valarray<double>(n);
		d.boundary_west  = valarray<double>(n);
		d.planNeumann();

		// solve over south interface, and save results
		valarray<double> north_block(n * n);
		valarray<double> east_block(n * n);
		valarray<double> south_block(n * n);
		valarray<double> west_block(n * n);
		for (int i = 0; i < n; i++) {
			d.boundary_south[i] = 1;
			d.solve();
			// fill the blocks

			north_block[slice(i * n, n, 1)] = d.u[slice(n * (n - 1), n, 1)];
			east_block[slice(i * n, n, 1)]  = d.u[slice((n - 1), n, n)];
			south_block[slice(i * n, n, 1)] = d.u[slice(0, n, 1)];
			west_block[slice(i * n, n, 1)]  = d.u[slice(0, n, n)];
			south_block[i * n + i] -= 1;
			d.boundary_south[i] = 0;
		}

		//now insert these results into the matrix for each interface
        for(Iface iface: todo){
			bool reverse_x = false;
			if ((iface.axis == X_AXIS && !iface.right) || (iface.axis == Y_AXIS && iface.right)) {
				reverse_x = true;
			}

			for (int x = 0; x < iface.l_south; x++) {
				int j       = iface.i_south + x;
				int block_i = x;
				if (reverse_x) {
					block_i = iface.l_south - 1 - x;
				}

				vector<double> row;
				vector<int>    global;
				row.reserve(n * 2 + n * 2);
				global.reserve(n * 3 + n * 4);
				for (int y = 0; y < iface.l_south; y++) {
					int block_j = y;
					if (reverse_x) {
						block_j = iface.l_south - 1 - y;
					}
					row.push_back(-south_block[block_i * n + block_j]);
					global.push_back(iface.i_south + y);
				}
				if (iface.i_west != -1) {
					for (int y = 0; y < iface.l_south; y++) {
						int block_j = y;
						row.push_back(-west_block[block_i * n + block_j]);
						if (iface.right) {
							global.push_back(iface.i_west + y);
						} else {
							global.push_back(iface.i_west + iface.l_south - 1 - y);
						}
					}
				}
				if (iface.i_north != -1) {
					for (int y = 0; y < iface.l_south; y++) {
						int block_j = y;
						if (reverse_x) {
							block_j = iface.l_south - 1 - y;
						}
						row.push_back(-north_block[block_i * n + block_j]);
						global.push_back(iface.i_north + y);
					}
				}
				if (iface.i_east != -1) {
					for (int y = 0; y < iface.l_south; y++) {
						int block_j = y;
						row.push_back(-east_block[block_i * n + block_j]);
						if (iface.right) {
							global.push_back(iface.i_east + y);
						} else {
							global.push_back(iface.i_east + iface.l_south - 1 - y);
						}
					}
				}
				// insert row for domain
				if (j == delete_row) {
					double one=1;
					A->insertGlobalValues(j, 1, &one, &j);
				} else {
					A->insertGlobalValues(j, row.size(), &row[0], &global[0]);
				}
			}
		}
	}

	// cerr << "Num types: " << num_types << "\n";
	// transpose matrix and return
	A->fillComplete();
    */
	return A;
}
RCP<matrix_type> DomainCollection::formInvDiag(RCP<map_type> map, int del_row)
{
	int              size = max(n, n);
	RCP<matrix_type> A    = rcp(new matrix_type(map, size * 6));
    /*
    // create iface objects
	set<pair<Iface,Iface>> ifaces;
	auto       iface_view = iface_info->getLocalView<Kokkos::HostSpace>();
	for (size_t i = 0; i < iface_view.dimension(0); i += 22) {
		Iface right;
		Iface left;

		right.right       = true;
		right.axis        = iface_view(i + 21, 0);
		right.global_i[0] = iface_view(i, 0);
		right.types[0]    = iface_view(i, 1);

		for (int q = 1; q < 4; q++) {
			right.global_i[q * 2] = iface_view(i + q * 3, 0);
			right.types[q]        = iface_view(i + q * 3 + 1, 0);
		}

		left.right   = false;
		left.axis    = iface_view(i + 21, 0);
		left.global_i[0] = iface_view(i, 0);
		left.types[0]    = iface_view(i, 1);
		for (int q = 1; q < 4; q++) {
			left.global_i[q * 2] = iface_view(i + 9 + q * 3, 0);
			left.types[q]        = iface_view(i + 9 + q * 3 + 1, 0);
		}

		pair<Iface,Iface> p;
        p.first = left;
        p.second = right;
		ifaces.insert(p);
	}

	while (!ifaces.empty()) {
		// the first in the set is the type of interface that we are going to solve for
		set<Iface> todo;
		pair<Iface, Iface> curr_pair = *ifaces.begin();

		ifaces.erase(ifaces.begin());
		Iface curr_type_left  = curr_pair.first;
		Iface curr_type_right = curr_pair.second;
		// todo.insert(curr_type_left);
		todo.insert(curr_type_right);
		bool contains_del_row = (del_row != -1 && curr_type_right.i_south <= del_row
		                         && del_row < curr_type_right.i_south + curr_type_right.l_south);
		if (!contains_del_row) {
			set<pair<Iface, Iface>> to_be_deleted;
			for (auto iter = ifaces.begin(); iter != ifaces.end(); iter++) {
				if (*iter == curr_pair) {
					todo.insert(iter->second);
					to_be_deleted.insert(*iter);
				}
			}
			for (auto i : to_be_deleted) {
				ifaces.erase(i);
			}
		}

		// create domain representing curr_type_left
		DomainSignature ds;
		if (curr_type_left.t_north == NEUMANN) {
			ds.nbr[0] = -1;
		} else {
			ds.nbr[0] = 1;
		}
		if (curr_type_left.t_east == NEUMANN) {
			ds.nbr[2] = -1;
		} else {
			ds.nbr[2] = 1;
		}
		if (curr_type_left.t_south == NEUMANN) {
			ds.nbr[4] = -1;
		} else {
			ds.nbr[4] = 1;
		}
		if (curr_type_left.t_west == NEUMANN) {
			ds.nbr[6] = -1;
		} else {
			ds.nbr[6] = 1;
		}
		Domain d_left(ds, n, n, h_x, h_y);
		d_left.boundary_north = valarray<double>(n);
		d_left.boundary_south = valarray<double>(n);
		d_left.boundary_east  = valarray<double>(n);
		d_left.boundary_west  = valarray<double>(n);
		d_left.planNeumann();

		if (curr_type_right.t_north == NEUMANN) {
			ds.nbr[0] = -1;
		} else {
			ds.nbr[0] = 1;
		}
		if (curr_type_right.t_east == NEUMANN) {
			ds.nbr[2] = -1;
		} else {
			ds.nbr[2] = 1;
		}
		if (curr_type_right.t_south == NEUMANN) {
			ds.nbr[4] = -1;
		} else {
			ds.nbr[4] = 1;
		}
		if (curr_type_right.t_west == NEUMANN) {
			ds.nbr[6] = -1;
		} else {
			ds.nbr[6] = 1;
		}
		Domain d_right(ds, n, n, h_x, h_y);
		d_right.boundary_north = valarray<double>(n);
		d_right.boundary_south = valarray<double>(n);
		d_right.boundary_east  = valarray<double>(n);
		d_right.boundary_west  = valarray<double>(n);
		d_right.planNeumann();

		// solve over south interface, and save results
		valarray<double> south_block(n * n);
		for (int i = 0; i < n; i++) {
			d_right.boundary_south[i] = 1;
            d_left.boundary_south[n-1-i] =1;
			d_right.solve();
			d_left.solve();

			// fill the blocks
			south_block[slice(i * n, n, 1)] = d_right.u[slice(0, n, 1)];
			for (int j = 0; j < n; j++) {
				south_block[i * n + j] += d_left.u[n - 1 - j];
			}
			south_block[i * n + i] -= 2;

			d_right.boundary_south[i] = 0;
            d_left.boundary_south[n-1-i] =0;
		}

		//south_block *= 2;
		int internal_i = del_row % n;
		if(contains_del_row){
            for(int x = 0;x<n;x++){
				south_block[internal_i * n + x] =0;
			}
			south_block[internal_i * n + internal_i] = -1;
		}

        //compute inverse of block
        valarray<int> ipiv(n+1);
        int lwork = n*n;
        valarray<double> work(lwork);
        int info;
		dgetrf_(&n, &n, &south_block[0], &n, &ipiv[0], &info);
		dgetri_(&n, &south_block[0], &n, &ipiv[0], &work[0], &lwork, &info);

		// now insert these results into the matrix for each interface
		for(Iface iface: todo){
			for (int x = 0; x < iface.l_south; x++) {
				int j = iface.i_south + x;
				vector<double> row;
				vector<int>    global;
				row.reserve(n * 2 + n * 2);
				global.reserve(n * 3 + n * 4);
				for (int y = 0; y < iface.l_south; y++) {
					row.push_back(-south_block[x * n + y]);
					global.push_back(iface.i_south + y);
				}
				// insert row for domain
				A->insertGlobalValues(j, row.size(), &row[0], &global[0]);
			}
		}
	}

	// transpose matrix and return
	A->fillComplete();
    */
	return A;
}
RCP<RBMatrix> DomainCollection::formRBMatrix(RCP<map_type> map, int delete_row)
{
	// create iface objects
	set<Iface> ifaces;
	auto       iface_view = iface_info->getLocalView<Kokkos::HostSpace>();
	for (size_t i = 0; i < iface_view.dimension(0); i += 22) {
		Iface right;
		Iface left;

		right.right       = true;
		right.axis        = iface_view(i + 21, 0);
		right.global_i[0] = iface_view(i, 0);
		right.types[0]    = iface_view(i+1, 1);

		for (int q = 1; q < 4; q++) {
			right.global_i[q * 2] = iface_view(i + q * 3, 0);
			right.types[q]        = iface_view(i + q * 3 + 1, 0);
		}

		left.right       = false;
		left.axis        = iface_view(i + 21, 0);
		left.global_i[0] = iface_view(i, 0);
		left.types[0]    = iface_view(i+1, 1);
		for (int q = 1; q < 4; q++) {
			left.global_i[q * 2] = iface_view(i + 9 + q * 3, 0);
			left.types[q]        = iface_view(i + 9 + q * 3 + 1, 0);
		}

		ifaces.insert(left);
		ifaces.insert(right);
	}

	int           size      = max(n, n);
	RCP<RBMatrix> A         = rcp(new RBMatrix(map, size, ifaces.size() / 2));
	int           num_types = 0;
	while (!ifaces.empty()) {
		num_types++;
		// the first in the set is the type of interface that we are going to solve for
		set<Iface> todo;
		Iface      curr_type = *ifaces.begin();
		ifaces.erase(ifaces.begin());
		todo.insert(curr_type);
		set<Iface> to_be_deleted;
		for (auto iter = ifaces.begin(); iter != ifaces.end(); iter++) {
			if (*iter == curr_type) {
				todo.insert(*iter);
				to_be_deleted.insert(*iter);

				// TODO fix this iterator
				// iter=ifaces.begin();
			}
		}
		for (Iface i : to_be_deleted) {
			ifaces.erase(i);
		}

		// create domain representing curr_type
		DomainSignature ds;
        for(int q=0;q<4;q++){
			if (curr_type.types[q] == NEUMANN) {
				ds.nbr[(q * 2 + 4) % 8] = -1;
			} else {
				ds.nbr[(q * 2 + 4) % 8] = 1;
			}
		}
		Domain d(ds, n, n, h_x, h_y);
		d.boundary_north = valarray<double>(n);
		d.boundary_south = valarray<double>(n);
		d.boundary_east  = valarray<double>(n);
		d.boundary_west  = valarray<double>(n);
		d.planNeumann();

		// solve over south interface, and save results
		RCP<valarray<double>> north_block_ptr = rcp(new valarray<double>(n * n));
		valarray<double> &    north_block     = *north_block_ptr;

		RCP<valarray<double>> east_block_ptr = rcp(new valarray<double>(n * n));
		valarray<double> &    east_block     = *east_block_ptr;

		RCP<valarray<double>> south_block_ptr = rcp(new valarray<double>(n * n));
		valarray<double> &    south_block     = *south_block_ptr;

		RCP<valarray<double>> west_block_ptr = rcp(new valarray<double>(n * n));
		valarray<double> &    west_block     = *west_block_ptr;

		for (int i = 0; i < n; i++) {
			d.boundary_south[i] = 1;
			d.solve();
			// fill the blocks

			north_block[slice(i*n, n, 1)] = d.getDiffNorth();
			east_block[slice(i*n, n, 1)]  = d.getDiffEast();
			south_block[slice(i*n, n, 1)] = d.getDiffSouth();
			west_block[slice(i*n, n, 1)]  = d.getDiffWest();
			d.boundary_south[i] = 0;
		}

		//now insert these results into the matrix for each interface
		for (Iface iface : todo) {
			bool reverse_x
			= (iface.axis == X_AXIS && !iface.right) || (iface.axis == Y_AXIS && iface.right);
			bool reverse_y = !iface.right;

			int j = iface.global_i[0];
			int i = iface.global_i[0];

			A->insertBlock(i, j, south_block_ptr, reverse_x, reverse_x);

			if (iface.global_i[2] != -1) {
				i = iface.global_i[2];
				A->insertBlock(i, j, west_block_ptr, reverse_y, reverse_x);
			}
			if (iface.global_i[4] != -1) {
				i = iface.global_i[4];
				A->insertBlock(i, j, north_block_ptr, reverse_x, reverse_x);
			}
			if (iface.global_i[6] != -1) {
				i = iface.global_i[6];
				A->insertBlock(i, j, east_block_ptr, reverse_y, reverse_x);
			}
		}
	}

	// cerr << "Num types: " << num_types << "\n";
	// transpose matrix and return
	//A->fillComplete();
	A->skip_index = delete_row;
	A->createRangeMap();
	return A;

}

void DomainCollection::outputSolution(std::ostream &os)
{
	int num_i, num_j, d_x;
	if (amr) {
		num_i = n * sqrt(num_global_domains / 5);
		num_j = n * sqrt(num_global_domains / 5);
		d_x   = sqrt(num_global_domains / 5);
	} else {
		num_i = n * sqrt(num_global_domains);
		num_j = n * sqrt(num_global_domains);
		d_x   = sqrt(num_global_domains);
	}
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = domain_i * d_x + domain_j;
			os << domains[id]->u[internal_i * n + internal_j] << '\n';
		}
	}
}
void DomainCollection::outputSolutionRefined(std::ostream &os)
{
	int num_i = 2 * n * sqrt(num_global_domains / 5);
	int num_j = 2 * n * sqrt(num_global_domains / 5);
	int d_x   = 2 * sqrt(num_global_domains / 5);
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = d_x * d_x / 4 + domain_i * d_x + domain_j;
			os << domains[id]->u[internal_i * n + internal_j] << '\n';
		}
	}
}
void DomainCollection::outputResidual(std::ostream &os)
{
	int num_i, num_j, d_x;
	if (amr) {
		num_i = n * sqrt(num_global_domains / 5);
		num_j = n * sqrt(num_global_domains / 5);
		d_x   = sqrt(num_global_domains / 5);
	} else {
		num_i = n * sqrt(num_global_domains);
		num_j = n * sqrt(num_global_domains);
		d_x   = sqrt(num_global_domains);
	}
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = domain_i * d_x + domain_j;
			os << domains[id]->resid[internal_i * n + internal_j] << '\n';
		}
	}
}
void DomainCollection::outputResidualRefined(std::ostream &os)
{
	int num_i = 2 * n * sqrt(num_global_domains / 5);
	int num_j = 2 * n * sqrt(num_global_domains / 5);
	int d_x   = 2 * sqrt(num_global_domains / 5);
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = d_x * d_x / 4 + domain_i * d_x + domain_j;
			os << domains[id]->resid[internal_i * n + internal_j] << '\n';
		}
	}
}
void DomainCollection::outputError(std::ostream &os)
{
	int num_i, num_j, d_x;
	if (amr) {
		num_i = n * sqrt(num_global_domains / 5);
		num_j = n * sqrt(num_global_domains / 5);
		d_x   = sqrt(num_global_domains / 5);
	} else {
		num_i = n * sqrt(num_global_domains);
		num_j = n * sqrt(num_global_domains);
		d_x   = sqrt(num_global_domains);
	}
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = domain_i * d_x + domain_j;
			os << domains[id]->error[internal_i * n + internal_j] << '\n';
		}
		}
}
void DomainCollection::outputErrorRefined(std::ostream &os)
{
	int num_i = 2 * n * sqrt(num_global_domains / 5);
	int num_j = 2 * n * sqrt(num_global_domains / 5);
	int d_x   = 2 * sqrt(num_global_domains / 5);
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = d_x * d_x / 4 + domain_i * d_x + domain_j;
			os << domains[id]->error[internal_i * n + internal_j] << '\n';
		}
	}
}
