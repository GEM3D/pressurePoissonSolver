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

DomainCollection::DomainCollection(DomainSignatureCollection dsc, int n, double h_x, double h_y,
                                   RCP<const Teuchos::Comm<int>> comm)
{
	// cerr<< "Low:  " << low << "\n";
	// cerr<< "High: " << high << "\n";
	this->comm         = comm;
	this->n            = n;
	this->h_x          = h_x;
	this->h_y          = h_y;
	num_global_domains = dsc.num_global_domains;
	for (auto p : dsc.domains) {
		DomainSignature ds       = p.second;
		int             i        = ds.id;

		// create a domain
		RCP<Domain> d_ptr = rcp(new Domain(ds, n, h_x, h_y));
		domains[i]        = d_ptr;
	}
}

void DomainCollection::initNeumann(function<double(double, double)> ffun,
                                   function<double(double, double)> efun,
                                   function<double(double, double)> nfunx,
                                   function<double(double, double)> nfuny, bool amr)
{
    this->amr = amr;
	for (auto &p : domains) {
		Domain &d        = *p.second;

		// Generate RHS vector
		std::valarray<double> &f     = d.f;
		std::valarray<double> &exact = d.exact;

		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x           = d.x_start + d.h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + d.h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = efun(x, y);
				// west
				if (!d.hasNbr(Side::west)) {
					d.boundary_west[yi] = nfunx(0.0, y);
				}
				// east
				if (!d.hasNbr(Side::east)) {
					double x = 1.0;
					if (amr) {
						x = 2.0;
					}
					d.boundary_east[yi] = nfunx(x, y);
				}
				// south
				if (!d.hasNbr(Side::south)) {
					d.boundary_south[xi] = nfuny(x, 0.0);
				}
				// north
				if (!d.hasNbr(Side::north)) {
					d.boundary_north[xi] = nfuny(x, 1.0);
				}
			}
		}
		d.planNeumann();
	}

	// create map for domains
	generateMaps();
	if (!amr) {
		distributeIfaceInfo();
	}
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
				double x           = d.x_start + d.h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + d.h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = gfun(x, y);
				if (!d.hasNbr(Side::north)) {
					d.boundary_north[xi] = gfun(x, 1.0);
				}
				if (!d.hasNbr(Side::east)) {
					d.boundary_east[yi] = gfun(2.0, y);
				}
				if (!d.hasNbr(Side::south)) {
					d.boundary_south[xi] = gfun(x, 0.0);
				}
				if (!d.hasNbr(Side::west)) {
					d.boundary_west[yi] = gfun(0.0, y);
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
				if (!d.hasNbr(Side::north)) {
					d.boundary_north[xi] = gfun(x, 1.0);
				}
				if (!d.hasNbr(Side::east)) {
					d.boundary_east[yi] = gfun(2.0, y);
				}
				if (!d.hasNbr(Side::south)) {
					d.boundary_south[xi] = gfun(x, 0.0);
				}
				if (!d.hasNbr(Side::west)) {
					d.boundary_west[yi] = gfun(0.0, y);
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
	int &       curr_i = num_cols;
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
			int              curr = queue.front();
			Domain &d    = *domains.at(curr);
			queue.pop_front();
			visited.insert(curr);
            not_visited.erase(curr);
            Side s = Side::north;
            do{
				if (d.hasNbr(s) && visited.count(d.nbr(s)) == 0) {
					// a new edge that we have not assigned an index to
					d.index(s) = curr_i;
					for (int i = 0; i < Iface::size; i++) {
						c_iface_global.push_back(curr_i*Iface::size + i);
						iface_global.push_back(curr_i*Iface::size + i);
						curr_c_i++;
					}
					for (int i = 0; i < n; i++) {
						global.push_back(curr_i*n + i);
						matrix_global.push_back(curr_i*n + i);
						curr_matrix_i++;
					}
					curr_i++;

					// fine case
					if (d.hasFineNbr(s)) {
						Domain &nbr_left  = *domains.at(d.nbr(s));
						Domain &nbr_right = *domains.at(d.nbrRight(s));

						// set center indexes
						nbr_left.indexCenter(!s)  = d.index(s);
						nbr_right.indexCenter(!s) = d.index(s);

						// set left and right indexes index
						nbr_left.index(!s) = curr_i;
						for (int i = 0; i < Iface::size; i++) {
							c_iface_global.push_back(curr_i * Iface::size + i);
							iface_global.push_back(curr_i * Iface::size + i);
							curr_c_i++;
						}
						for (int i = 0; i < n; i++) {
							global.push_back(curr_i * n + i);
							matrix_global.push_back(curr_i * n + i);
							curr_matrix_i++;
						}
						curr_i++;
						nbr_right.index(!s) = curr_i;
						for (int i = 0; i < Iface::size; i++) {
							c_iface_global.push_back(curr_i * Iface::size + i);
							iface_global.push_back(curr_i * Iface::size + i);
							curr_c_i++;
						}
						for (int i = 0; i < n; i++) {
							global.push_back(curr_i * n + i);
							matrix_global.push_back(curr_i * n + i);
							curr_matrix_i++;
						}
						curr_i++;

						// enqueue domains
						if (enqueued.count(d.nbr(s)) == 0) {
							queue.push_back(d.nbr(s));
							enqueued.insert(d.nbr(s));
						}
						if (enqueued.count(d.nbrRight(s)) == 0) {
							queue.push_back(d.nbrRight(s));
							enqueued.insert(d.nbrRight(s));
						}
						// coarse case
					} else if (d.hasCoarseNbr(s)) {
						// TODO
						// normal case
					} else {
						Domain &nbr = *domains.at(d.nbr(s));
						nbr.index(!s)        = d.index(s);
						// enqueue domain
						if (enqueued.count(d.nbr(s)) == 0) {
							queue.push_back(d.nbr(s));
							enqueued.insert(d.nbr(s));
						}
					}
				}
				s++;
			} while (s != Side::north);
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
		collection_map = Teuchos::rcp(new map_type(-1, &global[0], global.size(), 0, this->comm));
		collection_iface_map
		= Teuchos::rcp(new map_type(-1, &c_iface_global[0], c_iface_global.size(), 0, this->comm));
		matrix_map
		= Teuchos::rcp(new map_type(-1, &matrix_global[0], matrix_global.size(), 0, this->comm));
		iface_map
		= Teuchos::rcp(new map_type(-1, &iface_global[0], iface_global.size(), 0, this->comm));
	}
}
void DomainCollection::distributeIfaceInfo()
{
	int_vector_type dist(collection_iface_map, 1);
	iface_info = rcp(new int_vector_type(iface_map, 1));
	for (auto &p : domains) {
		Domain &d = *p.second;
		Iface::writeIfaces(d, dist);
	}
	Tpetra::Export<> exporter(collection_iface_map, iface_map);
	iface_info->doExport(dist, exporter, Tpetra::CombineMode::ADD);
	Iface::readIfaces(ifaces, *iface_info);
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
	diff.update(-2, gamma, 1);
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
	vector_type ghost(collection_map, 2);
	vector_type one_ghost(matrix_map, 2);
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
RCP<RBMatrix> DomainCollection::formRBMatrix(RCP<map_type> map, int delete_row)
{
	RCP<RBMatrix> A = rcp(new RBMatrix(map, n, num_cols));
	// create iface objects
	set<Iface> ifaces = this->ifaces;

	int num_types = 0;
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
				ds.nbr_id[(q * 2 + 4) % 8] = -1;
			} else {
				ds.nbr_id[(q * 2 + 4) % 8] = 1;
			}
		}
		Domain d(ds, n, h_x, h_y);
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
			north_block[slice(i * n, n, 1)] = d.boundary_north - d.getSide(Side::north);
			east_block[slice(i * n, n, 1)]  = d.boundary_east - d.getSide(Side::east);
			south_block[slice(i * n, n, 1)] = d.boundary_south - d.getSide(Side::south);
			west_block[slice(i * n, n, 1)]  = d.boundary_west - d.getSide(Side::west);

			d.boundary_south[i] = 0;
		}

		// now insert these results into the matrix for each interface
		for (Iface iface : todo) {
			bool reverse_x
			= (iface.axis == X_AXIS && iface.right) || (iface.axis == Y_AXIS && !iface.right);
			bool reverse_y = iface.right;

			int j = iface.global_i[0];
			int i = iface.global_i[0];

			A->insertBlock(i, j, south_block_ptr, reverse_x, reverse_x);

			if (iface.global_i[1] != -1) {
				i = iface.global_i[1];
				A->insertBlock(i, j, west_block_ptr, reverse_y, reverse_x);
			}
			if (iface.global_i[2] != -1) {
				i = iface.global_i[2];
				A->insertBlock(i, j, north_block_ptr, reverse_x, reverse_x);
			}
			if (iface.global_i[3] != -1) {
				i = iface.global_i[3];
				A->insertBlock(i, j, east_block_ptr, reverse_y, reverse_x);
			}
		}
	}

	// cerr << "Num types: " << num_types << "\n";
	// transpose matrix and return
	// A->fillComplete();
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
