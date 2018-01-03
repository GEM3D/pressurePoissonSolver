#include "SchurHelper.h"
#include <array>
#include <tuple>
using namespace std;
enum axis_enum { X_AXIS, Y_AXIS };
enum bc_enum { DIRICHLET, NEUMANN, REFINED };

SchurHelper::SchurHelper(DomainCollection dc, shared_ptr<PatchSolver> solver,
                         shared_ptr<PatchOperator> op, shared_ptr<Interpolator> interpolator)
{
	this->dc           = dc;
	this->solver       = solver;
	this->op           = op;
	this->interpolator = interpolator;
	local_gamma        = dc.getNewSchurDistVec();
	local_interp       = dc.getNewSchurDistVec();
	gamma              = dc.getNewSchurVec();
	PW<IS> dist_is;
	ISCreateBlock(MPI_COMM_WORLD, dc.n, dc.iface_dist_map_vec.size(), &dc.iface_dist_map_vec[0],
	              PETSC_COPY_VALUES, &dist_is);
	VecScatterCreate(gamma, dist_is, local_gamma, nullptr, &scatter);
}

void SchurHelper::solveWithInterface(const Vec f, Vec u, const Vec gamma, Vec diff)
{
	// initilize our local variables
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

    VecScale(local_interp,0);
	// solve over domains on this proc
	for (auto &p : dc.domains) {
		Domain &d = p.second;
		solver->solve(d, f, u, local_gamma);
		interpolator->interpolate(d, u, local_interp);
	}

	// export diff vector
	VecScale(diff, 0);
	VecScatterBegin(scatter, local_interp, diff, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, diff, ADD_VALUES, SCATTER_REVERSE);
	VecAXPBY(diff, 1.0, -1.0, gamma);
}
void SchurHelper::applyWithInterface(const Vec u, const Vec gamma, Vec f)
{
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
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
		set<MatrixBlock> row_blocks = iface.getGlobalRowBlocks();
		blocks.insert(row_blocks.begin(), row_blocks.end());
	}
	int num_types = 0;
	Vec u, f, r, e, gamma, interp;
	VecCreateSeq(PETSC_COMM_SELF, n * n, &u);
	VecCreateSeq(PETSC_COMM_SELF, n * n, &f);
	VecCreateSeq(PETSC_COMM_SELF, n * n, &r);
	VecCreateSeq(PETSC_COMM_SELF, n * n, &e);
	VecCreateSeq(PETSC_COMM_SELF, n, &gamma);
	VecCreateSeq(PETSC_COMM_SELF, n, &interp);
	double *interp_view, *gamma_view;
	VecGetArray(interp, &interp_view);
	VecGetArray(gamma, &gamma_view);
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
		ds.x_length        = curr_type.length;
		ds.y_length        = curr_type.length;
		ds.nbr_id[0]       = 1;
		ds.neumann         = curr_type.neumann;
		ds.zero_patch      = curr_type.zero_patch;
		ds.local_i         = {{0, 0, 0, 0}};
		ds.local_i_center  = {{0, 0, 0, 0}};
		ds.local_i_refined = {{0, 0, 0, 0, 0, 0, 0}};
		solver->addDomain(ds);

		map<BlockKey, shared_ptr<valarray<double>>> coeffs;
		// allocate blocks of coefficients
		for (const MatrixBlock &b : todo) {
			shared_ptr<valarray<double>> ptr = coeffs[b];
			if (ptr.get() == nullptr) {
				coeffs[b] = shared_ptr<valarray<double>>(new valarray<double>(n * n));
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
				VecScale(interp, 0);
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
	VecRestoreArray(interp, &interp_view);
	VecRestoreArray(gamma, &gamma_view);
	VecDestroy(&u);
	VecDestroy(&f);
	VecDestroy(&r);
	VecDestroy(&e);
	VecDestroy(&gamma);
	VecDestroy(&interp);
}
PW_explicit<Mat> SchurHelper::formCRSMatrix()
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
	PW<Mat> A;
	MatCreate(MPI_COMM_WORLD, &A);
	int local_size  = dc.ifaces.size() * n;
	int global_size = dc.num_global_interfaces * n;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 19 * n, nullptr, 19 * n, nullptr);

	auto insertBlock
	= [&](int i, int j, shared_ptr<valarray<double>> block, bool flip_i, bool flip_j) {
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
		  vector<int> inds_i(n);
		  iota(inds_i.begin(), inds_i.end(), local_i);
		  vector<int> inds_j(n);
		  iota(inds_j.begin(), inds_j.end(), local_j);

		  MatSetValues(A, n, &inds_i[0], n, &inds_j[0], &copy[0], ADD_VALUES);
	  };

	assembleMatrix(insertBlock);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	return A;
}
