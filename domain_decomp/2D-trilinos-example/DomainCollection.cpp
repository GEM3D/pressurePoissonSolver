#include "DomainCollection.h"
#include <tuple>
#include <array>
#include <Tpetra_Experimental_BlockCrsMatrix_def.hpp>
#include <Tpetra_Map_decl.hpp>
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

DomainCollection::DomainCollection(DomainSignatureCollection dsc, int n,
                                   RCP<const Teuchos::Comm<int>> comm)
{
	// cerr<< "Low:  " << low << "\n";
	// cerr<< "High: " << high << "\n";
	this->comm         = comm;
	this->n            = n;
    this->dsc = dsc;
	num_global_domains = dsc.num_global_domains;
	for (auto p : dsc.domains) {
		DomainSignature ds       = p.second;
		int             i        = ds.id;

		// create a domain
		RCP<Domain> d_ptr = rcp(new Domain(ds, n));
		domains[i]        = d_ptr;
	}
}

void DomainCollection::initNeumann(function<double(double, double)> ffun,
                                   function<double(double, double)> efun,
                                   function<double(double, double)> nfunx,
                                   function<double(double, double)> nfuny, bool amr)
{
	neumann = true;
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
					d.boundary_west[yi] = nfunx(d.x_start, y);
				}
				// east
				if (!d.hasNbr(Side::east)) {
					d.boundary_east[yi] = nfunx(d.x_start+d.x_length, y);
				}
				// south
				if (!d.hasNbr(Side::south)) {
					d.boundary_south[xi] = nfuny(x, d.y_start);
				}
				// north
				if (!d.hasNbr(Side::north)) {
					d.boundary_north[xi] = nfuny(x, d.y_start+d.y_length);
				}
			}
		}
		d.planNeumann();
	}

	// create map for domains
	generateMaps();
	distributeIfaceInfo();
	for (const Iface &i : ifaces) {
        const_cast<Iface&>(i).setNeumann();
	}
	getBlocks();
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
					d.boundary_north[xi] = gfun(x, d.y_start+d.y_length);
				}
				if (!d.hasNbr(Side::east)) {
					d.boundary_east[yi] = gfun(d.x_start+d.x_length, y);
				}
				if (!d.hasNbr(Side::south)) {
					d.boundary_south[xi] = gfun(x, d.y_start);
				}
				if (!d.hasNbr(Side::west)) {
					d.boundary_west[yi] = gfun(d.x_start, y);
				}
			}
		}
		d.planDirichlet();
	}
	// create map for domains
	generateMaps();
	distributeIfaceInfo();
    getBlocks();
}
void DomainCollection::getBlocks(){
	for (const Iface iface : ifaces) {
		bool reverse_x
		= (iface.axis == X_AXIS && iface.right) || (iface.axis == Y_AXIS && !iface.right);
		bool reverse_y = iface.right;

		int j = iface.global_i[0];
		if (iface.hasFineNbr[0]) {
			Blockk a(iface.global_i[0], j, reverse_x, reverse_x, iface.right, iface.neumann,
			         BlockType::north_coarse_in);
			if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
			Blockk b(iface.refined_left[0], j, reverse_x, reverse_x, iface.right, iface.neumann,
			         BlockType::north_coarse_out_left);
			if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
			Blockk c(iface.refined_right[0], j, reverse_x, reverse_x, iface.right, iface.neumann,
			         BlockType::north_coarse_out_right);
			if (matrix_map->getLocalElement(c.i * n) != -1) blocks.insert(c);
		} else if (iface.hasCoarseNbr[0]) {
			Blockk a(iface.global_i[0], j, reverse_x, reverse_x, iface.right, iface.neumann,
			         BlockType::north_fine_in);
			if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
			if (iface.isCoarseLeft[0]) {
				Blockk b(iface.center_i[0], j, reverse_x, reverse_x, iface.right, iface.neumann,
				         BlockType::north_fine_out_left);
				if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
			} else {
				Blockk b(iface.center_i[0], j, reverse_x, reverse_x, iface.right, iface.neumann,
				         BlockType::north_fine_out_right);
				if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
			}
		} else {
			Blockk a(iface.global_i[0], j, reverse_x, reverse_x, iface.right, iface.neumann,
			         BlockType::north);
			if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
		}

		if (iface.global_i[1] != -1) {
			if (iface.hasFineNbr[1]) {
				Blockk a(iface.global_i[1], j, reverse_y, reverse_x, iface.right, iface.neumann,
				         BlockType::east_coarse_in);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
				Blockk b(iface.refined_left[1], j, reverse_y, reverse_x, iface.right, iface.neumann,
				         BlockType::east_coarse_out_left);
				if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				Blockk c(iface.refined_right[1], j, reverse_y, reverse_x, iface.right,
				         iface.neumann, BlockType::east_coarse_out_right);
				if (matrix_map->getLocalElement(c.i * n) != -1) blocks.insert(c);
			} else if (iface.hasCoarseNbr[1]) {
				Blockk a(iface.global_i[1], j, reverse_y, reverse_x, iface.right, iface.neumann,
				         BlockType::east_fine_in);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
				if (iface.isCoarseLeft[1]) {
					Blockk b(iface.center_i[1], j, reverse_y, reverse_x, iface.right, iface.neumann,
					         BlockType::east_fine_out_left);
					if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				} else {
					Blockk b(iface.center_i[1], j, reverse_y, reverse_x, iface.right, iface.neumann,
					         BlockType::east_fine_out_right);
					if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				}
			} else {
				Blockk a(iface.global_i[1], j, reverse_y, reverse_x, iface.right, iface.neumann,
				         BlockType::east);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
			}
		}
		if (iface.global_i[2] != -1) {
			if (iface.hasFineNbr[2]) {
				Blockk a(iface.global_i[2], j, reverse_x, reverse_x, iface.right, iface.neumann,
				         BlockType::south_coarse_in);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
				Blockk b(iface.refined_left[2], j, reverse_x, reverse_x, iface.right, iface.neumann,
				         BlockType::south_coarse_out_left);
				if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				Blockk c(iface.refined_right[2], j, reverse_x, reverse_x, iface.right,
				         iface.neumann, BlockType::south_coarse_out_right);
				if (matrix_map->getLocalElement(c.i * n) != -1) blocks.insert(c);
			} else if (iface.hasCoarseNbr[2]) {
				Blockk a(iface.global_i[2], j, reverse_x, reverse_x, iface.right, iface.neumann,
				         BlockType::south_fine_in);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
				if (iface.isCoarseLeft[2]) {
					Blockk b(iface.center_i[2], j, reverse_x, reverse_x, iface.right, iface.neumann,
					         BlockType::south_fine_out_left);
					if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				} else {
					Blockk b(iface.center_i[2], j, reverse_x, reverse_x, iface.right, iface.neumann,
					         BlockType::south_fine_out_right);
					if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				}
			} else {
				Blockk a(iface.global_i[2], j, reverse_x, reverse_x, iface.right, iface.neumann,
				         BlockType::south);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
			}
		}
		if (iface.global_i[3] != -1) {
			if (iface.hasFineNbr[3]) {
				Blockk a(iface.global_i[3], j, reverse_y, reverse_x, iface.right, iface.neumann,
				         BlockType::west_coarse_in);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
				Blockk b(iface.refined_left[3], j, reverse_y, reverse_x, iface.right, iface.neumann,
				         BlockType::west_coarse_out_left);
				if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				Blockk c(iface.refined_right[3], j, reverse_y, reverse_x, iface.right,
				         iface.neumann, BlockType::west_coarse_out_right);
				if (matrix_map->getLocalElement(c.i * n) != -1) blocks.insert(c);
			} else if (iface.hasCoarseNbr[3]) {
				Blockk a(iface.global_i[3], j, reverse_y, reverse_x, iface.right, iface.neumann,
				         BlockType::west_fine_in);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
				if (iface.isCoarseLeft[3]) {
					Blockk b(iface.center_i[3], j, reverse_y, reverse_x, iface.right, iface.neumann,
					         BlockType::west_fine_out_left);
					if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				} else {
					Blockk b(iface.center_i[3], j, reverse_y, reverse_x, iface.right, iface.neumann,
					         BlockType::west_fine_out_right);
					if (matrix_map->getLocalElement(b.i * n) != -1) blocks.insert(b);
				}
			} else {
				Blockk a(iface.global_i[3], j, reverse_y, reverse_x, iface.right, iface.neumann,
				         BlockType::west);
				if (matrix_map->getLocalElement(a.i * n) != -1) blocks.insert(a);
			}
		}
	}
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
	int         curr_matrix_i = 0;
	global.reserve(domains.size() * (2 * n + 2 * n));
	auto addToMap = [&](int curr_i, int global_i) {
        if(global_i==-1){
            cerr << "neg global i"<<endl;
        }
		for (int i = 0; i < Iface::size; i++) {
			c_iface_global.push_back(global_i * Iface::size + i);
			curr_c_i++;
		}
		for (int i = 0; i < n; i++) {
			global.push_back(global_i * n + i);
			curr_matrix_i++;
		}
	};
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
				if (d.hasNbr(s) && d.index(s) == -1) {
					// a new edge that we have not assigned an index to
					d.index(s) = curr_i++;
					addToMap(d.index(s), d.globalIndex(s));

					// fine case
					if (d.hasFineNbr(s)) {
						d.indexRefinedLeft(s)  = curr_i++;
						d.indexRefinedRight(s) = curr_i++;
						addToMap(d.indexRefinedLeft(s), d.globalIndexRefinedLeft(s));
						addToMap(d.indexRefinedRight(s), d.globalIndexRefinedRight(s));

						// left
						try {
							Domain &nbr_left         = *domains.at(d.nbr(s));
							nbr_left.index(!s)       = d.indexRefinedLeft(s);
							nbr_left.indexCenter(!s) = d.index(s);
							if (enqueued.count(d.nbr(s)) == 0) {
								queue.push_back(d.nbr(s));
								enqueued.insert(d.nbr(s));
							}
						} catch (out_of_range &oor) {
							// do nothing
						}

						// right
						try {
							Domain &nbr_right         = *domains.at(d.nbrRight(s));
							nbr_right.index(!s)       = d.indexRefinedRight(s);
							nbr_right.indexCenter(!s) = d.index(s);
							if (enqueued.count(d.nbrRight(s)) == 0) {
								queue.push_back(d.nbrRight(s));
								enqueued.insert(d.nbrRight(s));
							}
						} catch (out_of_range &oor) {
							// do nothing
						}
						// coarse case
					} else if (d.hasCoarseNbr(s)) {
						d.indexCenter(s) = curr_i++;
						addToMap(d.indexCenter(s), d.globalIndexCenter(s));
                        int other_i = -1;
						try {
							Domain &nbr = *domains.at(d.nbr(s));
							other_i     = curr_i++;
							if (d.isCoarseLeft(s)) {
								nbr.indexRefinedLeft(!s)  = d.index(s);
								nbr.indexRefinedRight(!s) = other_i;
								addToMap(nbr.indexRefinedRight(!s), nbr.globalIndexRefinedRight(!s));
							} else {
								nbr.indexRefinedRight(!s) = d.index(s);
								nbr.indexRefinedLeft(!s)  = other_i;
								addToMap(nbr.indexRefinedLeft(!s), nbr.globalIndexRefinedLeft(!s));
							}
							nbr.index(!s) = d.indexCenter(s);

							// enqueue domains
							if (enqueued.count(nbr.ds.id) == 0) {
								queue.push_back(nbr.ds.id);
								enqueued.insert(nbr.ds.id);
							}
						} catch (out_of_range &oor) {
							// do nothing
						}
						try {
							int nbr_id = -1;
							if (d.isCoarseLeft(s)) {
								Side nbr_s = s;
								nbr_s--;
								nbr_id = d.nbr(nbr_s);
							} else {
								Side nbr_s = s;
								nbr_s++;
								nbr_id = d.nbr(nbr_s);
							}
							Domain &buddy        = *domains.at(nbr_id);
							buddy.indexCenter(s) = d.indexCenter(s);
							if (other_i == -1) {
								other_i = curr_i++;
								addToMap(other_i, buddy.globalIndex(s));
							}
							buddy.index(s) = other_i;
							// enqueue domains
							if (enqueued.count(buddy.ds.id) == 0) {
								queue.push_back(buddy.ds.id);
								enqueued.insert(buddy.ds.id);
							}
						} catch (out_of_range &oor) {
							// do nothing
						}
						// normal case
					} else {
						try {
							Domain &nbr   = *domains.at(d.nbr(s));
							nbr.index(!s) = d.index(s);
							// enqueue domain
							if (enqueued.count(d.nbr(s)) == 0) {
								queue.push_back(d.nbr(s));
								enqueued.insert(d.nbr(s));
							}
						} catch (out_of_range &oor) {
                            // do nothing
						}
					}
				}
				s++;
			} while (s != Side::north);
		}
	}

	vector<int> matrix_global;
	vector<int> iface_global;
	for (int i = n * dsc.matrix_j_low; i < n * dsc.matrix_j_high; i++) {
		matrix_global.push_back(i);
	}
	for (int i = 0; i < Iface::size * dsc.num_global_interfaces; i++) {
		iface_global.push_back(i);
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
		matrix_map = Teuchos::rcp(new map_type(-1, matrix_global.size(), 0, this->comm));
		iface_map = Tpetra::createLocalMap<int, int>(
		(size_t) Iface::size * dsc.num_global_interfaces, this->comm);
	}
#ifdef DNDEBUG
    auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
    matrix_map->describe(*out);
    collection_map->describe(*out);
#endif
}
void DomainCollection::distributeIfaceInfo()
{
	int_vector_type dist(collection_iface_map, 1);
	iface_info = rcp(new int_vector_type(iface_map, 1));
	for (auto &p : domains) {
		Domain &d = *p.second;
		Iface::writeIfaces(d, *iface_info);
	}
    iface_info->reduce();
	Iface::readIfaces(ifaces, *iface_info);
}
void DomainCollection::getFluxDiff(vector_type &diff)
{
	Tpetra::Import<> importer(diff.getMap(), collection_map);
	vector_type      local_diff(collection_map, 1);

	// solve over domains on this proc
	for (auto &p : domains) {
		p.second->getFluxDiff(local_diff);
	}

	diff.scale(0);
	diff.doExport(local_diff, importer, Tpetra::CombineMode::ADD);
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
	local_gamma.doImport(gamma, importer, Tpetra::CombineMode::INSERT);

	// solve over domains on this proc
	for (auto &p : domains) {
        Domain&d = *p.second;
		d.solveWithInterface(local_gamma);
	}

	if (neumann && zero_u) {
		//make avarage of solution be zero
		double avg = integrateU() / area();
		for (auto &p : domains) {
			Domain &d = *p.second;
			d.u -= avg;
		}
	}
	for (auto &p : domains) {
		p.second->getDiff(local_diff);
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
	// diff.update(-2, gamma, 1);
	diff.scale(-1);
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
double DomainCollection::integrateF()
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

	double sum = 0;
	for (auto &p : domains) {
		sum += p.second->integrateF();
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &sum, &retval);
	return retval;
}
double DomainCollection::integrateBoundaryFlux()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second->integrateBoundaryFlux();
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &sum, &retval);
	return retval;
}
double DomainCollection::area()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second->area();
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &sum, &retval);
	return retval;
}
double DomainCollection::integrateU()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second->integrateU();
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &sum, &retval);
	return retval;
}
double DomainCollection::integrateExact()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second->integrateExact();
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &sum, &retval);
	return retval;
}
double DomainCollection::integrateAU()
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

	double sum = 0;
	for (auto &p : domains) {
		sum += p.second->integrateAU();
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &sum, &retval);
	return retval;
}
void DomainCollection::formCrsMatrix(RCP<matrix_type> &A, RCP<single_vector_type> &s)
{
	set<int> rows;
	set<int> cols;
	for (Blockk b : blocks) {
		rows.insert(b.i);
		cols.insert(b.j);
	}

	rows.erase(-1);
	vector<int> rows_array;
	for (int i : rows) {
        cols.erase(i);
		for (int q = 0; q < n; q++)
			rows_array.push_back(i * n + q);
	}
	vector<int> cols_array(rows_array);
	for (int j : cols) {
		for (int q = 0; q < n; q++)
			cols_array.push_back(j * n + q);
	}
	RCP<map_type> row_map = rcp(new map_type(-1, &rows_array[0], rows_array.size(), 0, this->comm));
	RCP<map_type> col_map = rcp(new map_type(-1, &cols_array[0], cols_array.size(), 0, this->comm));

	A = rcp(new matrix_type(matrix_map, col_map, 5 * n));
	s = rcp(new single_vector_type(matrix_map));

	set<pair<int, int>> inserted;
	valarray<double> shift(n);
	valarray<double> shift_rev(n);
	auto insertBlock = [&](int i, int j, RCP<valarray<double>> block, bool flip_i, bool flip_j) {
		if (i == j) {
			if (flip_j) {
				for (int q = 0; q < n; q++) {
					s->sumIntoGlobalValue(j*n + q, shift_rev[q]);
				}
			} else {
				for (int q = 0; q < n; q++) {
					s->sumIntoGlobalValue(j*n + q, shift[q]);
				}
			}
		}
		int local_i = A->getRowMap()->getLocalElement(i * n);
		int local_j = A->getColMap()->getLocalElement(j * n);

		valarray<double> &orig = *block;
		valarray<double>  copy(n * n);
		for (int i = 0; i < n; i++) {
			int block_i = i;
			if (flip_i) {
				block_i = n - i - 1;
			}
			for (int j = 0; j < n; j++) {
				int block_j = j;
				if (flip_j) {
					block_j = n - j - 1;
				}
				copy[i * n + j] = orig[block_i * n + block_j];
			}
		}
		vector<int> inds(n);
		for (int q = 0; q < n; q++) {
			inds[q] = local_j + q;
		}
        if(inserted.count(make_pair(i,j))==0){
			for (int q = 0; q < n; q++) {
				A->insertLocalValues(local_i + q, n, &copy[q * n], &inds[0]);
			}
        }else{
            inserted.insert(make_pair(i,j));
			for (int q = 0; q < n; q++) {
				A->sumIntoLocalValues(local_i + q, n, &copy[q * n], &inds[0]);
			}
        }
	};

	// create iface objects
	set<Blockk> blocks = this->blocks;

	int num_types = 0;
	while (!blocks.empty()) {
		num_types++;
		// the first in the set is the type of interface that we are going to solve for
		set<Blockk> todo;
		Blockk      curr_type = *blocks.begin();
		blocks.erase(blocks.begin());
		todo.insert(curr_type);
		set<Blockk> to_be_deleted;
		for (auto iter = blocks.begin(); iter != blocks.end(); iter++) {
			if (*iter == curr_type) {
				todo.insert(*iter);
				to_be_deleted.insert(*iter);

				// TODO fix this iterator
				// iter=ifaces.begin();
			}
		}
		for (Blockk i : to_be_deleted) {
			blocks.erase(i);
		}

		// create domain representing curr_type
		DomainSignature ds;
        ds.x_length=n;
        ds.y_length=n;
        for(int q=0;q<4;q++){
			if (curr_type.neumann[q]) {
				ds.nbr_id[q * 2] = -1;
			} else {
				ds.nbr_id[q * 2] = 1;
			}
		}
		Domain d(ds, n);
		d.boundary_north = valarray<double>(n);
		d.boundary_east  = valarray<double>(n);
		d.boundary_south = valarray<double>(n);
		d.boundary_west  = valarray<double>(n);
		d.planNeumann();

		// solve over south interface, and save results
		RCP<valarray<double>> n_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> e_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> s_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> w_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    n_b   = *n_ptr;
		valarray<double> &    e_b   = *e_ptr;
		valarray<double> &    s_b   = *s_ptr;
		valarray<double> &    w_b   = *w_ptr;

        RCP<valarray<double>> nf_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ef_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sf_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wf_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    nf_b   = *nf_ptr;
		valarray<double> &    ef_b   = *ef_ptr;
		valarray<double> &    sf_b   = *sf_ptr;
		valarray<double> &    wf_b   = *wf_ptr;

        RCP<valarray<double>> nfl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> efl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sfl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wfl_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    nfl_b   = *nfl_ptr;
		valarray<double> &    efl_b   = *efl_ptr;
		valarray<double> &    sfl_b   = *sfl_ptr;
		valarray<double> &    wfl_b   = *wfl_ptr;

        RCP<valarray<double>> nfr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> efr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sfr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wfr_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    nfr_b   = *nfr_ptr;
		valarray<double> &    efr_b   = *efr_ptr;
		valarray<double> &    sfr_b   = *sfr_ptr;
		valarray<double> &    wfr_b   = *wfr_ptr;

        RCP<valarray<double>> nc_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ec_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sc_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wc_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    nc_b   = *nc_ptr;
		valarray<double> &    ec_b   = *ec_ptr;
		valarray<double> &    sc_b   = *sc_ptr;
		valarray<double> &    wc_b   = *wc_ptr;

        RCP<valarray<double>> ncl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ecl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> scl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wcl_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    ncl_b   = *ncl_ptr;
		valarray<double> &    ecl_b   = *ecl_ptr;
		valarray<double> &    scl_b   = *scl_ptr;
		valarray<double> &    wcl_b   = *wcl_ptr;

        RCP<valarray<double>> ncr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ecr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> scr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wcr_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    ncr_b   = *ncr_ptr;
		valarray<double> &    ecr_b   = *ecr_ptr;
		valarray<double> &    scr_b   = *scr_ptr;
		valarray<double> &    wcr_b   = *wcr_ptr;


		for (int i = 0; i < n; i++) {
			d.boundary_north[i] = 1;
			d.solve();

			shift[i]             = d.u.sum() * 2 / (num_global_domains * n * n);
			shift_rev[n - 1 - i] = shift[i];

			// fill the blocks
			n_b[slice(i, n, n)] = d.getDiff(Side::north);
			e_b[slice(i, n, n)] = d.getDiff(Side::east);
			s_b[slice(i, n, n)] = d.getDiff(Side::south);
			w_b[slice(i, n, n)] = d.getDiff(Side::west);

			nf_b[slice(i, n, n)] = d.getDiffFine(Side::north);
			ef_b[slice(i, n, n)] = d.getDiffFine(Side::east);
			sf_b[slice(i, n, n)] = d.getDiffFine(Side::south);
			wf_b[slice(i, n, n)] = d.getDiffFine(Side::west);

			nfl_b[slice(i, n, n)] = d.getDiffFineToCoarseLeft(Side::north);
			efl_b[slice(i, n, n)] = d.getDiffFineToCoarseLeft(Side::east);
			sfl_b[slice(i, n, n)] = d.getDiffFineToCoarseLeft(Side::south);
			wfl_b[slice(i, n, n)] = d.getDiffFineToCoarseLeft(Side::west);

			nfr_b[slice(i, n, n)] = d.getDiffFineToCoarseRight(Side::north);
			efr_b[slice(i, n, n)] = d.getDiffFineToCoarseRight(Side::east);
			sfr_b[slice(i, n, n)] = d.getDiffFineToCoarseRight(Side::south);
			wfr_b[slice(i, n, n)] = d.getDiffFineToCoarseRight(Side::west);

			nc_b[slice(i, n, n)] = d.getDiffCoarse(Side::north);
			ec_b[slice(i, n, n)] = d.getDiffCoarse(Side::east);
			sc_b[slice(i, n, n)] = d.getDiffCoarse(Side::south);
			wc_b[slice(i, n, n)] = d.getDiffCoarse(Side::west);

			ncl_b[slice(i, n, n)] = d.getDiffCoarseToFineLeft(Side::north);
			ecl_b[slice(i, n, n)] = d.getDiffCoarseToFineLeft(Side::east);
			scl_b[slice(i, n, n)] = d.getDiffCoarseToFineLeft(Side::south);
			wcl_b[slice(i, n, n)] = d.getDiffCoarseToFineLeft(Side::west);

			ncr_b[slice(i, n, n)] = d.getDiffCoarseToFineRight(Side::north);
			ecr_b[slice(i, n, n)] = d.getDiffCoarseToFineRight(Side::east);
			scr_b[slice(i, n, n)] = d.getDiffCoarseToFineRight(Side::south);
			wcr_b[slice(i, n, n)] = d.getDiffCoarseToFineRight(Side::west);

			d.boundary_north[i] = 0;
		}

		auto getBlock = [&](BlockType type) {
			RCP<valarray<double>> ret;
			switch (type) {
                case BlockType::north:
					ret = n_ptr;
					break;
                case BlockType::east:
					ret = e_ptr;
					break;
                case BlockType::south:
					ret = s_ptr;
					break;
                case BlockType::west:
					ret = w_ptr;
					break;
				// fine in
                case BlockType::north_fine_in:
					ret = nf_ptr;
					break;
                case BlockType::east_fine_in:
					ret = ef_ptr;
					break;
                case BlockType::south_fine_in:
					ret = sf_ptr;
					break;
                case BlockType::west_fine_in:
					ret = wf_ptr;
					break;
				// fine out
				// left
                case BlockType::north_fine_out_left:
					ret = nfl_ptr;
					break;
                case BlockType::east_fine_out_left:
					ret = efl_ptr;
					break;
                case BlockType::south_fine_out_left:
					ret = sfl_ptr;
					break;
                case BlockType::west_fine_out_left:
					ret = wfl_ptr;
					break;
				// right
                case BlockType::north_fine_out_right:
					ret = nfr_ptr;
					break;
                case BlockType::east_fine_out_right:
					ret = efr_ptr;
					break;
                case BlockType::south_fine_out_right:
					ret = sfr_ptr;
					break;
                case BlockType::west_fine_out_right:
					ret = wfr_ptr;
					break;
				// coarse in
                case BlockType::north_coarse_in:
					ret = nc_ptr;
					break;
                case BlockType::east_coarse_in:
					ret = ec_ptr;
					break;
                case BlockType::south_coarse_in:
					ret = sc_ptr;
					break;
                case BlockType::west_coarse_in:
					ret = wc_ptr;
					break;
				// coarse out
				// left
                case BlockType::north_coarse_out_left:
					ret = ncl_ptr;
					break;
                case BlockType::east_coarse_out_left:
					ret = ecl_ptr;
					break;
                case BlockType::south_coarse_out_left:
					ret = scl_ptr;
					break;
                case BlockType::west_coarse_out_left:
					ret = wcl_ptr;
					break;
				// right
                case BlockType::north_coarse_out_right:
					ret = ncr_ptr;
					break;
                case BlockType::east_coarse_out_right:
					ret = ecr_ptr;
					break;
                case BlockType::south_coarse_out_right:
					ret = scr_ptr;
					break;
                case BlockType::west_coarse_out_right:
					ret = wcr_ptr;
			}
			return ret;
		};
		// now insert these results into the matrix for each interface
		for (Blockk block : todo)
		{
			insertBlock(block.i, block.j, getBlock(block.type), block.flip_i, block.flip_j);
		}
	}

	A->fillComplete(matrix_map, matrix_map);
}
RCP<RBMatrix> DomainCollection::formRBMatrix(RCP<map_type> map, int delete_row)
{
	RCP<RBMatrix> A = rcp(new RBMatrix(map, n, dsc.matrix_j_high - dsc.matrix_j_low));
	A->skip_index=delete_row;
#if NDEBUG
	for (Iface i : ifaces) {
		cerr << i << endl;
	}
#endif
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
        ds.x_length=n;
        ds.y_length=n;
        for(int q=0;q<4;q++){
			if (curr_type.neumann[q]) {
				ds.nbr_id[q * 2] = -1;
			} else {
				ds.nbr_id[q * 2] = 1;
			}
		}
		Domain d(ds, n);
		d.boundary_north = valarray<double>(n);
		d.boundary_east  = valarray<double>(n);
		d.boundary_south = valarray<double>(n);
		d.boundary_west  = valarray<double>(n);
		d.planNeumann();

		// solve over south interface, and save results
		RCP<valarray<double>> n_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> e_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> s_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> w_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    n_b   = *n_ptr;
		valarray<double> &    e_b   = *e_ptr;
		valarray<double> &    s_b   = *s_ptr;
		valarray<double> &    w_b   = *w_ptr;

        RCP<valarray<double>> nf_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ef_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sf_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wf_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    nf_b   = *nf_ptr;
		valarray<double> &    ef_b   = *ef_ptr;
		valarray<double> &    sf_b   = *sf_ptr;
		valarray<double> &    wf_b   = *wf_ptr;

        RCP<valarray<double>> nfl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> efl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sfl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wfl_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    nfl_b   = *nfl_ptr;
		valarray<double> &    efl_b   = *efl_ptr;
		valarray<double> &    sfl_b   = *sfl_ptr;
		valarray<double> &    wfl_b   = *wfl_ptr;

        RCP<valarray<double>> nfr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> efr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sfr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wfr_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    nfr_b   = *nfr_ptr;
		valarray<double> &    efr_b   = *efr_ptr;
		valarray<double> &    sfr_b   = *sfr_ptr;
		valarray<double> &    wfr_b   = *wfr_ptr;

        RCP<valarray<double>> nc_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ec_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sc_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wc_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    nc_b   = *nc_ptr;
		valarray<double> &    ec_b   = *ec_ptr;
		valarray<double> &    sc_b   = *sc_ptr;
		valarray<double> &    wc_b   = *wc_ptr;

        RCP<valarray<double>> ncl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ecl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> scl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wcl_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    ncl_b   = *ncl_ptr;
		valarray<double> &    ecl_b   = *ecl_ptr;
		valarray<double> &    scl_b   = *scl_ptr;
		valarray<double> &    wcl_b   = *wcl_ptr;

        RCP<valarray<double>> ncr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ecr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> scr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wcr_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &    ncr_b   = *ncr_ptr;
		valarray<double> &    ecr_b   = *ecr_ptr;
		valarray<double> &    scr_b   = *scr_ptr;
		valarray<double> &    wcr_b   = *wcr_ptr;


        valarray<double> shift(n);
        valarray<double> shift_rev(n);
		for (int i = 0; i < n; i++) {
			d.boundary_north[i] = 1;
			d.solve();
			d.boundary_north[i] = 0;

			shift[i]             = d.u.sum() * 2 / (num_global_domains * n * n);
			shift_rev[n - 1 - i] = shift[i];

			// fill the blocks
			n_b[slice(i * n, n, 1)] = d.getDiff(Side::north);
			e_b[slice(i * n, n, 1)] = d.getDiff(Side::east);
			s_b[slice(i * n, n, 1)] = d.getDiff(Side::south);
			w_b[slice(i * n, n, 1)] = d.getDiff(Side::west);

			nf_b[slice(i * n, n, 1)] = d.getDiffFine(Side::north);
			ef_b[slice(i * n, n, 1)] = d.getDiffFine(Side::east);
			sf_b[slice(i * n, n, 1)] = d.getDiffFine(Side::south);
			wf_b[slice(i * n, n, 1)] = d.getDiffFine(Side::west);

			nfl_b[slice(i * n, n, 1)] = d.getDiffFineToCoarseLeft(Side::north);
			efl_b[slice(i * n, n, 1)] = d.getDiffFineToCoarseLeft(Side::east);
			sfl_b[slice(i * n, n, 1)] = d.getDiffFineToCoarseLeft(Side::south);
			wfl_b[slice(i * n, n, 1)] = d.getDiffFineToCoarseLeft(Side::west);

			nfr_b[slice(i * n, n, 1)] = d.getDiffFineToCoarseRight(Side::north);
			efr_b[slice(i * n, n, 1)] = d.getDiffFineToCoarseRight(Side::east);
			sfr_b[slice(i * n, n, 1)] = d.getDiffFineToCoarseRight(Side::south);
			wfr_b[slice(i * n, n, 1)] = d.getDiffFineToCoarseRight(Side::west);

			nc_b[slice(i * n, n, 1)] = d.getDiffCoarse(Side::north);
			ec_b[slice(i * n, n, 1)] = d.getDiffCoarse(Side::east);
			sc_b[slice(i * n, n, 1)] = d.getDiffCoarse(Side::south);
			wc_b[slice(i * n, n, 1)] = d.getDiffCoarse(Side::west);

			ncl_b[slice(i * n, n, 1)] = d.getDiffCoarseToFineLeft(Side::north);
			ecl_b[slice(i * n, n, 1)] = d.getDiffCoarseToFineLeft(Side::east);
			scl_b[slice(i * n, n, 1)] = d.getDiffCoarseToFineLeft(Side::south);
			wcl_b[slice(i * n, n, 1)] = d.getDiffCoarseToFineLeft(Side::west);

			ncr_b[slice(i * n, n, 1)] = d.getDiffCoarseToFineRight(Side::north);
			ecr_b[slice(i * n, n, 1)] = d.getDiffCoarseToFineRight(Side::east);
			scr_b[slice(i * n, n, 1)] = d.getDiffCoarseToFineRight(Side::south);
			wcr_b[slice(i * n, n, 1)] = d.getDiffCoarseToFineRight(Side::west);

		}

		// now insert these results into the matrix for each interface
		for (Iface iface : todo) {

			bool reverse_x
			= (iface.axis == X_AXIS && iface.right) || (iface.axis == Y_AXIS && !iface.right);
			bool reverse_y = iface.right;

			int j = iface.global_i[0];

			if (neumann) {
				if (reverse_x) {
					A->getShift(j) += shift_rev;
				} else {
					A->getShift(j) += shift;
				}
			}
			if (iface.hasFineNbr[0]) {
				A->insertBlock(iface.global_i[0], j, nc_ptr, reverse_x, reverse_x);
				A->insertBlock(iface.refined_left[0], j, ncl_ptr, reverse_x, reverse_x);
				A->insertBlock(iface.refined_right[0], j, ncr_ptr, reverse_x, reverse_x);
			} else if (iface.hasCoarseNbr[0]) {
				A->insertBlock(iface.global_i[0], j, nf_ptr, reverse_x, reverse_x);
				if (iface.isCoarseLeft[0]) {
					A->insertBlock(iface.center_i[0], j, nfl_ptr, reverse_x, reverse_x);
				} else {
					A->insertBlock(iface.center_i[0], j, nfr_ptr, reverse_x, reverse_x);
				}
			} else {
				A->insertBlock(iface.global_i[0], j, n_ptr, reverse_x, reverse_x);
			}

			if (iface.global_i[1] != -1) {
				if (iface.hasFineNbr[1]) {
					A->insertBlock(iface.global_i[1], j, ec_ptr, reverse_y, reverse_x);
					A->insertBlock(iface.refined_left[1], j, ecl_ptr, reverse_y, reverse_x);
					A->insertBlock(iface.refined_right[1], j, ecr_ptr, reverse_y, reverse_x);
				} else if (iface.hasCoarseNbr[1]) {
					A->insertBlock(iface.global_i[1], j, ef_ptr, reverse_y, reverse_x);
					if (iface.isCoarseLeft[1]) {
						A->insertBlock(iface.center_i[1], j, efl_ptr, reverse_y, reverse_x);
					} else {
						A->insertBlock(iface.center_i[1], j, efr_ptr, reverse_y, reverse_x);
					}
				} else {
					A->insertBlock(iface.global_i[1], j, e_ptr, reverse_y, reverse_x);
				}
			}
			if (iface.global_i[2] != -1) {
				if (iface.hasFineNbr[2]) {
					A->insertBlock(iface.global_i[2], j, sc_ptr, reverse_x, reverse_x);
					A->insertBlock(iface.refined_left[2], j, scl_ptr, reverse_x, reverse_x);
					A->insertBlock(iface.refined_right[2], j, scr_ptr, reverse_x, reverse_x);
				} else if (iface.hasCoarseNbr[2]) {
					A->insertBlock(iface.global_i[2], j, sf_ptr, reverse_x, reverse_x);
					if (iface.isCoarseLeft[2]) {
						A->insertBlock(iface.center_i[2], j, sfl_ptr, reverse_x, reverse_x);
					} else {
						A->insertBlock(iface.center_i[2], j, sfr_ptr, reverse_x, reverse_x);
					}
				} else {
					A->insertBlock(iface.global_i[2], j, s_ptr, reverse_x, reverse_x);
				}
			}
			if (iface.global_i[3] != -1) {
				if (iface.hasFineNbr[3]) {
					A->insertBlock(iface.global_i[3], j, wc_ptr, reverse_y, reverse_x);
					A->insertBlock(iface.refined_left[3], j, wcl_ptr, reverse_y, reverse_x);
					A->insertBlock(iface.refined_right[3], j, wcr_ptr, reverse_y, reverse_x);
				} else if (iface.hasCoarseNbr[3]) {
					A->insertBlock(iface.global_i[3], j, wf_ptr, reverse_y, reverse_x);
					if (iface.isCoarseLeft[3]) {
						A->insertBlock(iface.center_i[3], j, wfl_ptr, reverse_y, reverse_x);
					} else {
						A->insertBlock(iface.center_i[3], j, wfr_ptr, reverse_y, reverse_x);
					}
				} else {
					A->insertBlock(iface.global_i[3], j, w_ptr, reverse_y, reverse_x);
				}
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
            Domain &d = *domains[id];
			os << d.resid[internal_i * n + internal_j]*d.h_x*d.h_y << '\n';
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
            Domain &d = *domains[id];
			os << d.resid[internal_i * n + internal_j]*d.h_x*d.h_y << '\n';
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
void DomainCollection::outputClaw()
{
	ofstream     t_file("fort.t0000");
	const string tab = "\t";
	t_file << 0.0 << tab << "time" << endl;
	t_file << 3 << tab << "meqn" << endl;
	t_file << num_global_domains << tab << "ngrids" << endl;
	t_file << 2 << tab << "num_aux" << endl;
	t_file << 2 << tab << "num_dim" << endl;
	t_file.close();
	ofstream q_file("fort.q0000");

	q_file.precision(10);
	q_file << scientific;
	for (auto &p : domains) {
		Domain &d = *p.second;
		d.outputClaw(q_file);
	}
	q_file.close();
}
