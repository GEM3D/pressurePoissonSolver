#include "DomainCollection.h"
#include <Tpetra_Experimental_BlockCrsMatrix_def.hpp>
#include <array>
#include <tuple>
using Teuchos::RCP;
using Teuchos::rcp;
using namespace std;
enum axis_enum { X_AXIS, Y_AXIS };
enum bc_enum { DIRICHLET, NEUMANN, REFINED };

void DomainCollection::setPatchSolver(RCP<PatchSolver> psolver) { solver = psolver; }
DomainCollection::DomainCollection(DomainSignatureCollection     dsc,
                                   RCP<const Teuchos::Comm<int>> comm)
{
	this->comm         = comm;
	this->dsc          = dsc;
}

void DomainCollection::solveWithInterface(const vector_type &f, vector_type &u,
                                          const vector_type &gamma, vector_type &diff)
{
	// initilize our local variables
	diff.update(1, gamma, 0);
	Tpetra::Import<> importer(diff.getMap(), dsc.getSchurDistMap());
	Tpetra::Export<> exporter(dsc.getSchurDistMap(), diff.getMap());
	vector_type      local_gamma(dsc.getSchurDistMap(), 1);
	vector_type      local_diff(dsc.getSchurDistMap(), 1);
	vector_type      local_interp(dsc.getSchurDistMap(), 1);
	local_gamma.doImport(gamma, importer, Tpetra::CombineMode::INSERT);

	// solve over domains on this proc
	for (auto &p : dsc.domains) {
		DomainSignature &d = p.second;
		solver->solve(d, f, u, local_gamma);
		interpolator->interpolate(d, u, local_interp);
	}

	// export diff vector
	diff.scale(0);
	diff.doExport(local_interp, importer, Tpetra::CombineMode::ADD);
	diff.update(2.0, gamma, -2.0);
}
void DomainCollection::applyWithInterface(const vector_type &u, const vector_type &gamma,
                                          vector_type &f)
{
	Tpetra::Import<> importer(gamma.getMap(), dsc.getSchurDistMap());
	vector_type      local_gamma(dsc.getSchurDistMap(), 1);
	local_gamma.putScalar(0);
	local_gamma.doImport(gamma, importer, Tpetra::CombineMode::INSERT);
	for (auto &p : dsc.domains) {
		DomainSignature &d = p.second;
		op->apply(d, u, local_gamma, f);
	}
}
double DomainCollection::integrate(const vector_type &u)
{
	double sum    = 0;
	auto   u_view = u.get1dView();

	for (auto &p : dsc.domains) {
		DomainSignature &d     = p.second;
		int              start = d.n * d.n * d.id;

		double patch_sum = 0;
		for (int i = 0; i < d.n * d.n; i++) {
			patch_sum += u_view[start + i];
		}

		double h_x = d.x_length / d.n;
		double h_y = d.y_length / d.n;
		patch_sum *= h_x * h_y;

		sum += patch_sum;
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &sum, &retval);
	return retval;
}
double DomainCollection::area()
{
	double sum = 0;
	for (auto &p : dsc.domains) {
		DomainSignature &d = p.second;
		sum += d.x_length * d.y_length;
	}
	double retval;
	Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM, 1, &sum, &retval);
	return retval;
}
void DomainCollection::assembleMatrix(inserter insertBlock)
{
	int              n = dsc.n;
	set<MatrixBlock> blocks;
	for (auto p : dsc.ifaces) {
		Iface &          iface      = p.second;
		set<MatrixBlock> row_blocks = iface.getRowBlocks();
		blocks.insert(row_blocks.begin(), row_blocks.end());
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
void DomainCollection::formCRSMatrix(Teuchos::RCP<map_type> map, Teuchos::RCP<matrix_type> &A)
{
	int         n = dsc.n;
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

	A = rcp(new matrix_type(dsc.getSchurRowMap(), col_map, 5 * n));

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

	assembleMatrix(insertBlock);
	A->fillComplete(map, map);
}
