#include "DomainCollection.h"
using Teuchos::RCP;
using Teuchos::rcp;
using namespace std;
DomainCollection::DomainCollection(int low, int high, int nx, int ny, int d_x, int d_y, double h_x,
                                   double h_y, RCP<const Teuchos::Comm<int>> comm,
                                   function<double(double, double)>          ffun,
                                   function<double(double, double)>          gfun)
{
	// cerr<< "Low:  " << low << "\n";
	// cerr<< "High: " << high << "\n";
	this->comm         = comm;
	this->nx           = nx;
	this->ny           = ny;
	this->h_x          = h_x;
	this->h_y          = h_y;
	num_domains        = high - low;
	num_global_domains = d_x * d_y;
	domains            = map<int, Domain *>();
	for (int i = low; i <= high; i++) {
		int domain_x = i % d_y;
		int domain_y = i / d_x;

		// Generate RHS vector
		std::valarray<double> f(nx * ny);
		std::valarray<double> exact(nx * ny);
		// Use local indices to access the entries of f_data.
		for (int yi = 0; yi < ny; yi++) {
			for (int xi = 0; xi < nx; xi++) {
				int    index_x      = domain_x * nx + xi;
				int    index_y      = domain_y * ny + yi;
				double x            = h_x / 2.0 + 1.0 * index_x / (nx * d_x);
				double y            = h_y / 2.0 + 1.0 * index_y / (ny * d_y);
				f[yi * nx + xi]     = ffun(x, y);
				exact[yi * nx + xi] = gfun(x, y);
				if (index_x == 0) {
					f[yi * nx + xi] += -2.0 / (h_x * h_x) * gfun(0.0, y);
				}
				if (index_x == d_x * nx - 1) {
					f[yi * nx + xi] += -2.0 / (h_x * h_x) * gfun(1.0, y);
				}
				if (index_y == 0) {
					f[yi * nx + xi] += -2.0 / (h_y * h_y) * gfun(x, 0.0);
				}
				if (index_y == d_y * ny - 1) {
					f[yi * nx + xi] += -2.0 / (h_y * h_y) * gfun(x, 1.0);
				}
			}
		}
		// create a domain
		Domain *d_ptr = new Domain(f, exact, nx, ny, h_x, h_y);
		Domain &d     = *d_ptr;

		// determine its neighbors
		// north
		int ns_start_i = d_y * (d_x - 1) * ny;
		if (domain_y != d_y - 1) {
			int nbr_y        = domain_y + 1;
			d.nbr_north      = nbr_y * d_x + domain_x;
			d.global_i_north = (domain_y * d_x + domain_x) * nx + ns_start_i;
		}
		// east
		if (domain_x != d_x - 1) {
			int nbr_x       = domain_x + 1;
			d.nbr_east      = domain_y * d_x + nbr_x;
			d.global_i_east = (domain_x * d_x + domain_y) * ny;
		}
		// south
		if (domain_y != 0) {
			int nbr_y        = domain_y - 1;
			d.nbr_south      = nbr_y * d_x + domain_x;
			d.global_i_south = ((domain_y - 1) * d_x + domain_x) * nx + ns_start_i;
		}
		// west
		if (domain_x != 0) {
			int nbr_x       = domain_x - 1;
			d.nbr_west      = domain_y * d_x + nbr_x;
			d.global_i_west = ((domain_x - 1) * d_x + domain_y) * ny;
		}
		domains[i] = d_ptr;
	}

	// create map for domains
    generateMaps();
}
void DomainCollection::generateMaps()
{
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	int        first = domains.begin()->first;
	queue.push_back(first);
	enqueued.insert(first);
	int         curr_i = 0;
	int         curr_matrix_i = 0;
	vector<int> global;
    vector<int> matrix_global;
	global.reserve(num_domains * (2 * nx + 2 * ny));
	while (!queue.empty()) {
		int     curr = queue.front();
		Domain &d    = *domains[curr];
		queue.pop_front();
		visited.insert(curr);
		if (d.nbr_north != -1 && visited.count(d.nbr_north) == 0) {
			// a new edge that we have not assigned an index to
			try {
				Domain &nbr       = *domains.at(d.nbr_north);
				nbr.local_i_south = curr_i;
				if (enqueued.count(d.nbr_north) == 0) {
					queue.push_back(d.nbr_north);
					enqueued.insert(d.nbr_north);
				}
			} catch (const out_of_range &oor) {
				// do nothing
			}
			d.local_i_north = curr_i;
			for (int i = 0; i < nx; i++) {
				global.push_back(d.global_i_north + i);
				matrix_global.push_back(d.global_i_north + i);
                curr_matrix_i++;
			}
			curr_i += nx;
		}
		if (d.nbr_east != -1 && visited.count(d.nbr_east) == 0) {
			// a new edge that we have not assigned an index to
			try {
				Domain &nbr      = *domains.at(d.nbr_east);
				nbr.local_i_west = curr_i;
				if (enqueued.count(d.nbr_east) == 0) {
					queue.push_back(d.nbr_east);
					enqueued.insert(d.nbr_east);
				}
			} catch (const out_of_range &oor) {
				// do nothing
			}
			d.local_i_east = curr_i;
			for (int i = 0; i < ny; i++) {
				global.push_back(d.global_i_east + i);
				matrix_global.push_back(d.global_i_east + i);
                curr_matrix_i++;
			}
			curr_i += ny;
		}
		if (d.nbr_south != -1 && visited.count(d.nbr_south) == 0) {
			// a new edge that we have not assigned an index to
			try {
				Domain &nbr       = *domains.at(d.nbr_south);
				nbr.local_i_north = curr_i;
				if (enqueued.count(d.nbr_south) == 0) {
					queue.push_back(d.nbr_south);
					enqueued.insert(d.nbr_south);
				}
			} catch (const out_of_range &oor) {
				// do nothing
			}
			d.local_i_south = curr_i;
			for (int i = 0; i < nx; i++) {
				global.push_back(d.global_i_south + i);
			}
			curr_i += nx;
		}
		if (d.nbr_west != -1 && visited.count(d.nbr_west) == 0) {
			// a new edge that we have not assigned an index to
			try {
				Domain &nbr      = *domains.at(d.nbr_west);
				nbr.local_i_east = curr_i;
				if (enqueued.count(d.nbr_west) == 0) {
					queue.push_back(d.nbr_west);
					enqueued.insert(d.nbr_west);
				}
			} catch (const out_of_range &oor) {
				// do nothing
			}
			d.local_i_west = curr_i;
			for (int i = 0; i < ny; i++) {
				global.push_back(d.global_i_west + i);
			}
			curr_i += ny;
		}
	}
	// Now that the global indices have been calculated, we can create a map for the interface
	// points
	if (num_global_domains == 1) {
		// this is a special case for when there is only one domain
		collection_map = Teuchos::rcp(new map_type(1, 0, comm));
		matrix_map     = Teuchos::rcp(new map_type(1, 0, comm));
	} else {
		collection_map = Teuchos::rcp(new map_type(-1, &global[0], curr_i, 0, this->comm));
		matrix_map
		= Teuchos::rcp(new map_type(-1, &matrix_global[0], curr_matrix_i, 0, this->comm));
	}
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
double DomainCollection::exactNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second->exactNorm(), 2);
	}
	return sqrt(result);
}
RCP<matrix_type> DomainCollection::formMatrix(RCP<map_type> map)
{
	// create domain for forming matrix
	std::valarray<double> f(nx * ny);
	Domain                d(f, f, nx, ny, h_x, h_y);
	d.nbr_north           = 1;
	d.nbr_east            = 1;
	d.nbr_south           = 1;
	d.nbr_west            = 1;
	d.boundary_north      = valarray<double>(nx);
	d.boundary_south      = valarray<double>(nx);
	d.boundary_east       = valarray<double>(ny);
	d.boundary_west       = valarray<double>(ny);
	int              size = max(nx, ny);
	RCP<matrix_type> A    = rcp(new matrix_type(map, size * 6));
	// north boundary
	for (int i = 0; i < nx; i++) {
		d.boundary_north[i] = 1;
		d.solve();
		// create row and insert for each domain
		for (auto &p : domains) {
			Domain &d2 = *p.second;
			if (d2.nbr_north != -1) {
				vector<double> row;
				vector<int>    global;
				row.reserve(nx * 2 + ny * 2);
				global.reserve(nx * 3 + ny * 4);
				int global_row = d2.global_i_north + i;
				for (int i = 0; i < nx; i++) {
					row.push_back(-d.u[nx * (ny - 1) + i]);
					global.push_back(d2.global_i_north + i);
				}
				row[i] += 1;
				if (d2.nbr_east != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(-d.u[(i + 1) * nx - 1]);
						global.push_back(d2.global_i_east + i);
					}
				}
				if (d2.nbr_south != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(-d.u[i]);
						global.push_back(d2.global_i_south + i);
					}
				}
				if (d2.nbr_west != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(-d.u[i * nx]);
						global.push_back(d2.global_i_west + i);
					}
				}
				// insert row for domain
				A->insertGlobalValues(global_row, row.size(), &row[0], &global[0]);
			}
		}
		d.boundary_north[i] = 0;
	}
	// east boundary
	for (int i = 0; i < ny; i++) {
		d.boundary_east[i] = 1;
		d.solve();
		// create row and insert for each domain
		for (auto &p : domains) {
			Domain &d2 = *p.second;
			if (d2.nbr_east != -1) {
				vector<double> row;
				vector<int>    global;
				row.reserve(nx * 2 + ny * 2);
				global.reserve(nx * 3 + ny * 4);
				int global_row = d2.global_i_east + i;
				for (int i = 0; i < ny; i++) {
					row.push_back(-d.u[(i + 1) * nx - 1]);
					global.push_back(d2.global_i_east + i);
				}
				row[i] += 1;
				if (d2.nbr_north != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(-d.u[nx * (ny - 1) + i]);
						global.push_back(d2.global_i_north + i);
					}
				}
				if (d2.nbr_south != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(-d.u[i]);
						global.push_back(d2.global_i_south + i);
					}
				}
				if (d2.nbr_west != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(-d.u[i * nx]);
						global.push_back(d2.global_i_west + i);
					}
				}
				// insert row for domain
				A->insertGlobalValues(global_row, row.size(), &row[0], &global[0]);
			}
		}
		d.boundary_east[i] = 0;
	}
	// south boundary
	for (int i = 0; i < nx; i++) {
		d.boundary_south[i] = 1;
		d.solve();
		// create row and insert for each domain
		for (auto &p : domains) {
			Domain &d2 = *p.second;
			if (d2.nbr_south != -1) {
				vector<double> row;
				vector<int>    global;
				row.reserve(nx * 2 + ny * 2);
				global.reserve(nx * 3 + ny * 4);
				int global_row = d2.global_i_south + i;
				for (int i = 0; i < nx; i++) {
					row.push_back(-d.u[i]);
					global.push_back(d2.global_i_south + i);
				}
				row[i] += 1;
				if (d2.nbr_north != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(-d.u[nx * (ny - 1) + i]);
						global.push_back(d2.global_i_north + i);
					}
				}
				if (d2.nbr_east != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(-d.u[(i + 1) * nx - 1]);
						global.push_back(d2.global_i_east + i);
					}
				}
				if (d2.nbr_west != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(-d.u[i * nx]);
						global.push_back(d2.global_i_west + i);
					}
				}
				// insert row for domain
				A->insertGlobalValues(global_row, row.size(), &row[0], &global[0]);
			}
		}
		d.boundary_south[i] = 0;
	}
	// west boundary
	for (int i = 0; i < ny; i++) {
		d.boundary_west[i] = 1;
		d.solve();
		// create row and insert for each domain
		for (auto &p : domains) {
			Domain &d2 = *p.second;
			if (d2.nbr_west != -1) {
				vector<double> row;
				vector<int>    global;
				row.reserve(nx * 2 + ny * 2);
				global.reserve(nx * 3 + ny * 4);
				int global_row = d2.global_i_west + i;
				for (int i = 0; i < ny; i++) {
					row.push_back(-d.u[i * nx]);
					global.push_back(d2.global_i_west + i);
				}
				row[i] += 1;
				if (d2.nbr_north != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(-d.u[nx * (ny - 1) + i]);
						global.push_back(d2.global_i_north + i);
					}
				}
				if (d2.nbr_east != -1) {
					for (int i = 0; i < ny; i++) {
						row.push_back(-d.u[(i + 1) * nx - 1]);
						global.push_back(d2.global_i_east + i);
					}
				}
				if (d2.nbr_south != -1) {
					for (int i = 0; i < nx; i++) {
						row.push_back(-d.u[i]);
						global.push_back(d2.global_i_south + i);
					}
				}
				// insert row for domain
				A->insertGlobalValues(global_row, row.size(), &row[0], &global[0]);
			}
		}
		d.boundary_west[i] = 0;
	}
	// transpose matrix and return
	A->fillComplete();
	return A;
}
