#include "SchurHelper.h"
#include <array>
#include <tuple>
using Teuchos::RCP;
using Teuchos::rcp;
using namespace std;
enum axis_enum { X_AXIS, Y_AXIS };
enum bc_enum { DIRICHLET, NEUMANN, REFINED };

SchurHelper::SchurHelper(DomainCollection dc, Teuchos::RCP<const Teuchos::Comm<int>> comm,
                         Teuchos::RCP<PatchSolver> solver, Teuchos::RCP<PatchOperator> op,
                         Teuchos::RCP<Interpolator> interpolator)
{
	this->dc           = dc;
	this->comm         = comm;
	this->solver       = solver;
	this->op           = op;
	this->interpolator = interpolator;
}

void SchurHelper::solveWithInterface(const vector_type &f, vector_type &u, const vector_type &gamma,
                                     vector_type &diff)
{
	// initilize our local variables
	diff.update(1, gamma, 0);
	Tpetra::Import<> importer(diff.getMap(), dc.getSchurDistMap());
	vector_type      local_gamma(dc.getSchurDistMap(), 1);
	vector_type      local_interp(dc.getSchurDistMap(), 1);
	local_gamma.doImport(gamma, importer, Tpetra::CombineMode::INSERT);

	// solve over domains on this proc
	for (auto &p : dc.domains) {
		Domain &d = p.second;
		solver->solve(d, f, u, local_gamma);
		interpolator->interpolate(d, u, local_interp);
	}

	// export diff vector
	diff.scale(0);
	diff.doExport(local_interp, importer, Tpetra::CombineMode::ADD);
	diff.update(1.0, gamma, -1.0);
}
void SchurHelper::applyWithInterface(const vector_type &u, const vector_type &gamma, vector_type &f)
{
	Tpetra::Import<> importer(gamma.getMap(), dc.getSchurDistMap());
	vector_type      local_gamma(dc.getSchurDistMap(), 1);
	local_gamma.putScalar(0);
	local_gamma.doImport(gamma, importer, Tpetra::CombineMode::INSERT);
	for (auto &p : dc.domains) {
		Domain &d = p.second;
		op->apply(d, u, local_gamma, f);
	}
}

struct BlockKey {
	InterpCase type;
	Side       s;

	BlockKey() {}
	BlockKey(const MatrixBlock &b)
	{
		s    = b.s;
		type = b.type;
	}
	friend bool operator<(const BlockKey &l, const BlockKey &r)
	{
		return std::tie(l.s, l.type) < std::tie(r.s, r.type);
	}
};
void SchurHelper::assembleMatrix(inserter insertBlock)
{
	int              n = dc.n;
	set<MatrixBlock> blocks;
	for (auto p : dc.ifaces) {
		Iface &          iface      = p.second;
		set<MatrixBlock> row_blocks = iface.getRowBlocks();
		blocks.insert(row_blocks.begin(), row_blocks.end());
	}
	int                 num_types       = 0;
	RCP<const map_type> local_map       = Tpetra::createLocalMap<int, int>(n * n, comm);
	RCP<const map_type> gamma_local_map = Tpetra::createLocalMap<int, int>(n, comm);
	vector_type         u(local_map, 1);
	vector_type         f(local_map, 1);
	vector_type         gamma(gamma_local_map, 1);
	vector_type         interp(gamma_local_map, 1);
	auto                interp_view = interp.get1dView();
	auto                gamma_view  = gamma.get1dViewNonConst();
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
			}
		}
		for (MatrixBlock i : to_be_deleted) {
			blocks.erase(i);
		}

		// create domain representing curr_type
		Domain ds;
		ds.n               = n;
		ds.id              = 0;
		ds.id_local        = 0;
		ds.x_length        = n;
		ds.y_length        = n;
		ds.nbr_id[0]       = 1;
		ds.neumann         = curr_type.neumann;
		ds.zero_patch      = curr_type.zero_patch;
		ds.local_i         = {{0, 0, 0, 0}};
		ds.local_i_center  = {{0, 0, 0, 0}};
		ds.local_i_refined = {{0, 0, 0, 0, 0, 0, 0}};
		solver->addDomain(ds);

		map<BlockKey, RCP<valarray<double>>> coeffs;
		// allocate blocks of coefficients
		for (const MatrixBlock &b : todo) {
			RCP<valarray<double>> ptr = coeffs[b];
			if (ptr.is_null()) {
				coeffs[b] = rcp(new valarray<double>(n * n));
			}
		}

		for (int j = 0; j < n; j++) {
			gamma_view[j] = 1;
			solver->solve(ds, f, u, gamma);
			gamma_view[j] = 0;

			// fill the blocks
			for (auto &p : coeffs) {
				Side       s    = p.first.s;
				InterpCase type = p.first.type;
				interp.putScalar(0);
				interpolator->interpolate(ds, s, type, u, interp);
				valarray<double> &block = *p.second;
				for (int i = 0; i < n; i++) {
					block[i * n + j] = -interp_view[i];
				}
				if (s == Side::north) {
					switch (type) {
						case InterpCase::normal:
							block[n * j + j] += 0.5;
							break;
						case InterpCase::coarse_from_coarse:
						case InterpCase::fine_from_fine_on_left:
						case InterpCase::fine_from_fine_on_right:
							block[n * j + j] += 1;
							break;
						default:
							break;
					}
				}
			}
		}

		// now insert these results into the matrix for each interface
		for (MatrixBlock block : todo) {
			insertBlock(block.i, block.j, coeffs[block], block.flip_i, block.flip_j);
		}
	}
}
Teuchos::RCP<matrix_type> SchurHelper::formCRSMatrix()
{
	int         n = dc.n;
	vector<int> cols_array;
	for (int i : dc.iface_map_vec) {
		for (int j = 0; j < n; j++) {
			cols_array.push_back(i * n + j);
		}
	}
	for (int i : dc.iface_off_proc_map_vec) {
		for (int j = 0; j < n; j++) {
			cols_array.push_back(i * n + j);
		}
	}
	RCP<map_type> col_map = rcp(new map_type(-1, &cols_array[0], cols_array.size(), 0, this->comm));

	RCP<matrix_type> A = rcp(new matrix_type(dc.getSchurRowMap(), col_map, 5 * n));

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
	A->fillComplete(dc.getSchurRowMap(), dc.getSchurRowMap());
	return A;
}
