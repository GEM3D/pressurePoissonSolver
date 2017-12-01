#include "DomainCollection.h"
#include <Tpetra_Experimental_BlockCrsMatrix_def.hpp>
#include <array>
#include <tuple>
#ifdef HAVE_VTK
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLPMultiBlockDataWriter.h>
#endif
using Teuchos::RCP;
using Teuchos::rcp;
using namespace std;
enum axis_enum { X_AXIS, Y_AXIS };
enum bc_enum { DIRICHLET, NEUMANN, REFINED };

void DomainCollection::setPatchSolver(RCP<PatchSolver> psolver){
    solver=psolver;
}
DomainCollection::DomainCollection(DomainSignatureCollection dsc, int n,
                                   RCP<const Teuchos::Comm<int>> comm)
{
	this->comm         = comm;
	this->n            = n;
	this->dsc          = dsc;
	num_global_domains = dsc.num_global_domains;
	for (auto p : dsc.domains) {
		DomainSignature ds = p.second;
		int             i  = ds.id;

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
	neumann   = true;
	this->amr = amr;
	for (auto &p : domains) {
		Domain &d = *p.second;

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
					d.boundary_east[yi] = nfunx(d.x_start + d.x_length, y);
				}
				// south
				if (!d.hasNbr(Side::south)) {
					d.boundary_south[xi] = nfuny(x, d.y_start);
				}
				// north
				if (!d.hasNbr(Side::north)) {
					d.boundary_north[xi] = nfuny(x, d.y_start + d.y_length);
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
					d.boundary_north[xi] = gfun(x, d.y_start + d.y_length);
				}
				if (!d.hasNbr(Side::east)) {
					d.boundary_east[yi] = gfun(d.x_start + d.x_length, y);
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
}

void DomainCollection::generateMaps()
{
}
void DomainCollection::distributeIfaceInfo() {}
RCP<vector_type> DomainCollection::getInterfaceCoords()
{
	RCP<vector_type> xy = rcp(new vector_type(dsc.getSchurRowMap(n), 2));
	Tpetra::Import<> importer(dsc.getSchurRowMap(n), dsc.getSchurDistMap(n));
	vector_type      local_xy(dsc.getSchurDistMap(n), 2);

	// solve over domains on this proc
	for (auto &p : domains) {
		p.second->putCoords(local_xy);
	}

	xy->doExport(local_xy, importer, Tpetra::CombineMode::ADD);
	return xy;
}
void DomainCollection::solveWithInterface(const vector_type &f, vector_type &u,const vector_type &gamma, vector_type &diff)
{
    // initilize our local variables
	diff.update(1, gamma, 0);
	Tpetra::Import<> importer(diff.getMap(), dsc.getSchurDistMap(n));
	Tpetra::Export<> exporter(dsc.getSchurDistMap(n), diff.getMap());
	vector_type      local_gamma(dsc.getSchurDistMap(n), 1);
	vector_type      local_diff(dsc.getSchurDistMap(n), 1);
	vector_type      local_interp(dsc.getSchurDistMap(n), 1);
	local_gamma.doImport(gamma, importer, Tpetra::CombineMode::INSERT);

	// solve over domains on this proc
	for (auto &p : domains) {
		Domain &d = *p.second;
        solver->solve(d.ds,f,u,local_gamma);
		d.solveWithInterface(local_gamma);
	}
    
    // interplate to interface points
	for (auto &p : domains) {
        interpolator->interpolate(p.second->ds,u,local_interp);
	}

	// export diff vector
	diff.scale(0);
	diff.doExport(local_interp, importer, Tpetra::CombineMode::ADD);
	diff.update(2.0,gamma,-2.0);
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
	vector_type ghost(dsc.getSchurDistMap(n), 2);
	vector_type one_ghost(dsc.getSchurRowMap(n), 2);
	for (auto &p : domains) {
		p.second->putGhostCells(ghost);
	}
	Tpetra::Export<> exporter(dsc.getSchurDistMap(n), dsc.getSchurRowMap(n));
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
	vector_type ghost(dsc.getSchurDistMap(n), 2);
	vector_type one_ghost(dsc.getSchurRowMap(n), 2);
	for (auto &p : domains) {
		p.second->putGhostCells(ghost);
	}
	Tpetra::Export<> exporter(dsc.getSchurDistMap(n), dsc.getSchurRowMap(n));
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
	vector_type ghost(dsc.getSchurDistMap(n), 2);
	vector_type one_ghost(dsc.getSchurRowMap(n), 2);
	for (auto &p : domains) {
		p.second->putGhostCells(ghost);
	}
	Tpetra::Export<> exporter(dsc.getSchurDistMap(n), dsc.getSchurRowMap(n));
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
void DomainCollection::assembleMatrix(inserter insertBlock,
                                     int n)
{
	set<MatrixBlock> blocks;
	for (auto p : dsc.ifaces) {
		Iface &          iface      = p.second;
		set<MatrixBlock> row_blocks = iface.getRowBlocks();
		blocks.insert(row_blocks.begin(), row_blocks.end());
	}
	if (n == -1) {
		n = this->n;
	}
	int num_types = 0;
	while (!blocks.empty()) {
		num_types++;
		// the first in the set is the type of interface that we are going to solve for
		set<MatrixBlock> todo;
		MatrixBlock      curr_type = *blocks.begin();
		blocks.erase(blocks.begin());
		todo.insert(curr_type);
		set<MatrixBlock> to_be_deleted;
		for (auto iter = blocks.begin(); iter != blocks.end(); iter++) {
			if (*iter == curr_type) {
				todo.insert(*iter);
				to_be_deleted.insert(*iter);

				// TODO fix this iterator
				// iter=ifaces.begin();
			}
		}
		for (MatrixBlock i : to_be_deleted) {
			blocks.erase(i);
		}

		// create domain representing curr_type
		DomainSignature ds;
		ds.x_length = n;
		ds.y_length = n;
		for (int q = 0; q < 4; q++) {
			if (curr_type.neumann[q]) {
				ds.nbr_id[q * 2] = -1;
			} else {
				ds.nbr_id[q * 2] = 1;
			}
		}
		ds.neumann    = curr_type.neumann;
		ds.zero_patch = curr_type.zero_patch;
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

		valarray<double> &n_b = *n_ptr;
		valarray<double> &e_b = *e_ptr;
		valarray<double> &s_b = *s_ptr;
		valarray<double> &w_b = *w_ptr;

		RCP<valarray<double>> nf_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ef_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sf_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wf_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &nf_b = *nf_ptr;
		valarray<double> &ef_b = *ef_ptr;
		valarray<double> &sf_b = *sf_ptr;
		valarray<double> &wf_b = *wf_ptr;

		RCP<valarray<double>> nfl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> efl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sfl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wfl_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &nfl_b = *nfl_ptr;
		valarray<double> &efl_b = *efl_ptr;
		valarray<double> &sfl_b = *sfl_ptr;
		valarray<double> &wfl_b = *wfl_ptr;

		RCP<valarray<double>> nfr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> efr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sfr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wfr_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &nfr_b = *nfr_ptr;
		valarray<double> &efr_b = *efr_ptr;
		valarray<double> &sfr_b = *sfr_ptr;
		valarray<double> &wfr_b = *wfr_ptr;

		RCP<valarray<double>> nc_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ec_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> sc_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wc_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &nc_b = *nc_ptr;
		valarray<double> &ec_b = *ec_ptr;
		valarray<double> &sc_b = *sc_ptr;
		valarray<double> &wc_b = *wc_ptr;

		RCP<valarray<double>> ncl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ecl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> scl_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wcl_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &ncl_b = *ncl_ptr;
		valarray<double> &ecl_b = *ecl_ptr;
		valarray<double> &scl_b = *scl_ptr;
		valarray<double> &wcl_b = *wcl_ptr;

		RCP<valarray<double>> ncr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> ecr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> scr_ptr = rcp(new valarray<double>(n * n));
		RCP<valarray<double>> wcr_ptr = rcp(new valarray<double>(n * n));

		valarray<double> &ncr_b = *ncr_ptr;
		valarray<double> &ecr_b = *ecr_ptr;
		valarray<double> &scr_b = *scr_ptr;
		valarray<double> &wcr_b = *wcr_ptr;

		for (int i = 0; i < n; i++) {
			d.boundary_north[i] = 1;
			d.solve();

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

		auto getBlock = [&](const MatrixBlock &b) {
			RCP<valarray<double>> ret;
			switch (b.type) {
				case BlockType::plain:
					switch (b.s) {
						case Side::north:
							ret = n_ptr;
							break;
						case Side::east:
							ret = e_ptr;
							break;
						case Side::south:
							ret = s_ptr;
							break;
						case Side::west:
							ret = w_ptr;
					}
					break;
				// fine in
				case BlockType::fine:
					switch (b.s) {
						case Side::north:
							ret = nf_ptr;
							break;
						case Side::east:
							ret = ef_ptr;
							break;
						case Side::south:
							ret = sf_ptr;
							break;
						case Side::west:
							ret = wf_ptr;
					}
					break;
				// fine out
				// left
				case BlockType::fine_out_left:
					switch (b.s) {
						case Side::north:
							ret = nfl_ptr;
							break;
						case Side::east:
							ret = efl_ptr;
							break;
						case Side::south:
							ret = sfl_ptr;
							break;
						case Side::west:
							ret = wfl_ptr;
					}
					break;
				// right
				case BlockType::fine_out_right:
					switch (b.s) {
						case Side::north:
							ret = nfr_ptr;
							break;
						case Side::east:
							ret = efr_ptr;
							break;
						case Side::south:
							ret = sfr_ptr;
							break;
						case Side::west:
							ret = wfr_ptr;
					}
					break;

				// coarse in
				case BlockType::coarse:
					switch (b.s) {
						case Side::north:
							ret = nc_ptr;
							break;
						case Side::east:
							ret = ec_ptr;
							break;
						case Side::south:
							ret = sc_ptr;
							break;
						case Side::west:
							ret = wc_ptr;
					}
					break;
				// coarse out
				// left
				case BlockType::coarse_out_left:
					switch (b.s) {
						case Side::north:
							ret = ncl_ptr;
							break;
						case Side::east:
							ret = ecl_ptr;
							break;
						case Side::south:
							ret = scl_ptr;
							break;
						case Side::west:
							ret = wcl_ptr;
					}
					break;
				// right
				case BlockType::coarse_out_right:
					switch (b.s) {
						case Side::north:
							ret = ncr_ptr;
							break;
						case Side::east:
							ret = ecr_ptr;
							break;
						case Side::south:
							ret = scr_ptr;
							break;
						case Side::west:
							ret = wcr_ptr;
					}
			}
			return ret;
		};
		// now insert these results into the matrix for each interface
		for (MatrixBlock block : todo) {
			insertBlock(block.i, block.j, getBlock(block), block.flip_i, block.flip_j);
		}
	}

}
void DomainCollection::formCRSMatrix(Teuchos::RCP<map_type> map, Teuchos::RCP<matrix_type> &A,
                                     int n)
{
	if (n == -1) {
		n = this->n;
	}
	vector<int> cols_array;
	for (int i : dsc.iface_map_vec) {
		for (int j = 0; j < n; j++) {
			cols_array.push_back(i * n + j);
		}
	}
	for (int i : dsc.iface_off_proc_map_vec) {
		for (int j = 0; j < n; j++) {
			cols_array.push_back(i * n + j);
		}
	}
	RCP<map_type> col_map = rcp(new map_type(-1, &cols_array[0], cols_array.size(), 0, this->comm));

	A  = rcp(new matrix_type(dsc.getSchurRowMap(n), col_map, 5 * n));

	set<pair<int, int>> inserted;
	auto insertBlock = [&](int i, int j, RCP<valarray<double>> block, bool flip_i, bool flip_j) {
		int local_i = i * n;
		int local_j = j * n;

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
		if (inserted.count(make_pair(i, j)) == 0) {
			for (int q = 0; q < n; q++) {
				A->insertLocalValues(local_i + q, n, &copy[q * n], &inds[0]);
			}
		} else {
			inserted.insert(make_pair(i, j));
			for (int q = 0; q < n; q++) {
				A->sumIntoLocalValues(local_i + q, n, &copy[q * n], &inds[0]);
			}
		}
	};

    assembleMatrix(insertBlock,n);
	A->fillComplete(map, map);
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
			int     domain_i   = i / n;
			int     internal_i = i % n;
			int     id         = domain_i * d_x + domain_j;
			Domain &d          = *domains[id];
			os << d.resid[internal_i * n + internal_j] * d.h_x * d.h_y << '\n';
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
			int     domain_i   = i / n;
			int     internal_i = i % n;
			int     id         = d_x * d_x / 4 + domain_i * d_x + domain_j;
			Domain &d          = *domains[id];
			os << d.resid[internal_i * n + internal_j] * d.h_x * d.h_y << '\n';
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
#ifdef HAVE_VTK
void DomainCollection::outputVTK()
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	set<vtkSmartPointer<vtkImageData>>   images;
	set<vtkSmartPointer<vtkDoubleArray>> arrays;
	// create MultiPieceDataSet and fill with patch information
	vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer
	= vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
	vtkSmartPointer<vtkMultiBlockDataSet> block = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	vtkSmartPointer<vtkMultiPieceDataSet> data  = vtkSmartPointer<vtkMultiPieceDataSet>::New();

	data->SetNumberOfPieces(domains.size());

	int i = 0;
	for (auto &p : domains) {
		Domain &d = *p.second;

		// create image object
		vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
		images.insert(image);
		image->SetOrigin(d.x_start, d.y_start, 0);
		image->SetSpacing(d.h_x, d.h_y, 0);
		image->SetExtent(d.x_start, d.x_start + d.x_length, d.y_start, d.y_start + d.y_length, 0,
		                 0);
		image->PrepareForNewData();
		image->SetDimensions(n + 1, n + 1, 1);

		// create solution vector
		vtkSmartPointer<vtkDoubleArray> solution = vtkSmartPointer<vtkDoubleArray>::New();
		arrays.insert(solution);
		solution->SetNumberOfComponents(1);
		solution->SetNumberOfValues(d.u.size());
		solution->SetName("Solution");
		for (size_t i = 0; i < d.u.size(); i++) {
			solution->SetValue(i, d.u[i]);
		}
		image->GetCellData()->AddArray(solution);

		// create error vector
		vtkSmartPointer<vtkDoubleArray> error = vtkSmartPointer<vtkDoubleArray>::New();
		arrays.insert(error);
		error->SetNumberOfComponents(1);
		error->SetNumberOfValues(d.error.size());
		error->SetName("Error");
		for (size_t i = 0; i < d.error.size(); i++) {
			error->SetValue(i, d.error[i]);
		}
		image->GetCellData()->AddArray(error);

		// create residual vector
		vtkSmartPointer<vtkDoubleArray> resid = vtkSmartPointer<vtkDoubleArray>::New();
		arrays.insert(resid);
		resid->SetNumberOfComponents(1);
		resid->SetNumberOfValues(d.resid.size());
		resid->SetName("Residual");
		for (size_t i = 0; i < d.resid.size(); i++) {
			resid->SetValue(i, d.resid[i] * d.h_x * d.h_y);
		}
		image->GetCellData()->AddArray(resid);

		// add image to dataset
		data->SetPiece(i, image);
		i++;
	}

	block->SetNumberOfBlocks(1);
	block->SetBlock(0, data);

	writer->SetFileName("blah.vtmb");
	// writer->SetNumberOfPieces(1);
	// writer->SetStartPiece(0);
	writer->SetInputData(block);
	writer->Update();

	// write data
	writer->Write();
}
#endif
