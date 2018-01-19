#include "SchurHelper.h"
#include "PBMatrix.h"
#include <array>
#include <iostream>
#include <numeric>
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
	ISCreateBlock(MPI_COMM_WORLD, dc.n * dc.n, dc.iface_dist_map_vec.size(),
	              &dc.iface_dist_map_vec[0], PETSC_COPY_VALUES, &dist_is);
	VecScatterCreate(gamma, dist_is, local_gamma, nullptr, &scatter);
}

void SchurHelper::solveWithInterface(const Vec f, Vec u, const Vec gamma, Vec diff)
{
	// initilize our local variables
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	VecScale(local_interp, 0);
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
void SchurHelper::solveWithSolution(const Vec f, Vec u)
{
	// initilize our local variables
	VecScale(local_gamma, 0);
	/*
	for (auto &p : dc.domains) {
	    Domain &d = p.second;
	    interpolator->interpolate(d, u, local_gamma);
	}
	*/
	/*
	VecScale(gamma, 0);
	VecScatterBegin(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	*/

	// solve over domains on this proc
	for (auto &p : dc.domains) {
		Domain &d = p.second;
		solver->solve(d, f, u, local_gamma);
	}
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
void SchurHelper::apply(const Vec u, Vec f)
{
	VecScale(local_interp, 0);
	for (auto &p : dc.domains) {
		Domain &d = p.second;
		interpolator->interpolate(d, u, local_interp);
	}
	VecScale(gamma, 0);
	VecScatterBegin(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	for (auto &p : dc.domains) {
		Domain &d = p.second;
		op->apply(d, u, local_gamma, f);
	}
}
enum class Rotation : char { x_cw, x_ccw, y_cw, y_ccw, z_cw, z_ccw };
constexpr bool sideIsLeft(const Side s)
{
	return (s == Side::north || s == Side::west || s == Side::bottom);
}
struct Block {
	static const Side             side_table[6][6];
	static const char             rots_table[6][6];
	static const vector<Rotation> main_rot_plan[6];
	char                          data = 0;
	Side                          main;
	Side                          aux;
	int                           j;
	int                           i;
	bitset<6>                     neumann;
	Block(Side main, int j, Side aux, int i, bitset<6> neumann)
	{
		this->main    = main;
		this->j       = j;
		this->aux     = aux;
		this->i       = i;
		this->neumann = neumann;
		data          = 0;
		bool left     = sideIsLeft(main);
		data |= left << 7;
		left = sideIsLeft(aux);
		data |= left << 3;
		rotate();
	}
	Block &operator*=(const Rotation &rot)
	{
		char r;
		// main rotation
		r = (data & ~(~0u << 2) << 4) >> 4;
		r = (r + rots_table[static_cast<int>(rot)][static_cast<int>(main)]) & 0b11;
		data &= ~(~(~0u << 2) << 4);
		data |= r << 4;
		// aux rotation
		r = data & ~(~0u << 2);
		r = (r + rots_table[static_cast<int>(rot)][static_cast<int>(aux)]) & 0b11;
		data &= ~0u << 2;
		data |= r;
		main                  = side_table[static_cast<int>(rot)][static_cast<int>(main)];
		aux                   = side_table[static_cast<int>(rot)][static_cast<int>(aux)];
		bitset<6> old_neumann = neumann;
		for (int i = 0; i < 6; i++) {
			neumann[(int) side_table[(int) rot][i]] = old_neumann[i];
		}
		return *this;
	}
	void rotate()
	{
		for (Rotation rot : main_rot_plan[static_cast<int>(main)]) {
			(*this) *= rot;
		}
		cerr << "NEUMANN: " << neumann << endl;
		cerr << "AUX:     " << (int) aux << endl;
	}
	bool operator<(const Block &b) const
	{
		return std::tie(i, j, data) < std::tie(b.i, b.j, b.data);
	}
	bool operator==(const Block &b) const { return neumann.to_ulong()==b.neumann.to_ulong(); }
	//
	bool mainLeft() { return sideIsLeft(main); }
	bool mainFlipped()
	{
		bool left      = sideIsLeft(main);
		bool orig_left = (data >> 7) & 0b1;
		return left != orig_left;
	}
	int  mainRot() { return (data >> 4) & 0b11; }
	bool auxLeft() { return sideIsLeft(aux); }
	bool auxFlipped()
	{
		bool left      = sideIsLeft(aux);
		bool orig_left = (data >> 3) & 0b1;
		return left != orig_left;
	}
	int auxRot() { return data & 0b11; }
};
struct BlockKey {
	Side          s;
	unsigned char neumann;

	BlockKey() {}
	BlockKey(const Block &b)
	{
		s       = b.aux;
		neumann = b.neumann.to_ulong();
	}
	friend bool operator<(const BlockKey &l, const BlockKey &r)
	{
		return l.s < r.s;
	}
};

const Side Block::side_table[6][6]
= {{Side::west, Side::east, Side::top, Side::bottom, Side::south, Side::north},
   {Side::west, Side::east, Side::bottom, Side::top, Side::north, Side::south},
   {Side::bottom, Side::top, Side::south, Side::north, Side::east, Side::west},
   {Side::top, Side::bottom, Side::south, Side::north, Side::west, Side::east},
   {Side::north, Side::south, Side::west, Side::east, Side::bottom, Side::top},
   {Side::south, Side::north, Side::east, Side::west, Side::bottom, Side::top}};
const char Block::rots_table[6][6] = {{3, 1, 0, 0, 2, 2}, {1, 3, 2, 2, 0, 0}, {1, 3, 3, 1, 1, 3},
                                      {1, 3, 1, 3, 3, 1}, {0, 0, 0, 0, 3, 1}, {0, 0, 0, 0, 1, 3}};
const vector<Rotation> Block::main_rot_plan[6] = {{},
                                                  {Rotation::z_cw, Rotation::z_cw},
                                                  {Rotation::z_cw},
                                                  {Rotation::z_ccw},
                                                  {Rotation::y_ccw},
                                                  {Rotation::y_cw}};
void SchurHelper::assembleMatrix(inserter insertBlock)
{
	int        n = dc.n;
	set<Block> blocks;
	for (auto &p : dc.ifaces) {
		IfaceSet &ifs = p.second;
		int       i   = ifs.id_global;
		for (const Iface &iface : ifs.ifaces) {
			Side aux = iface.s;
			for (int s = 0; s < 6; s++) {
				int j = iface.global_id[s];
				if (j != -1) {
					Side main = static_cast<Side>(s);
					blocks.insert(Block(main, j, aux, i, iface.neumann));
				}
			}
		}
	}
	int num_types = 0;
	Vec u, f, r, e, gamma, interp;
	VecCreateSeq(PETSC_COMM_SELF, n * n * n, &u);
	VecCreateSeq(PETSC_COMM_SELF, n * n * n, &f);
	VecCreateSeq(PETSC_COMM_SELF, n * n * n, &r);
	VecCreateSeq(PETSC_COMM_SELF, n * n * n, &e);
	VecCreateSeq(PETSC_COMM_SELF, n * n, &gamma);
	VecCreateSeq(PETSC_COMM_SELF, n * n, &interp);
	double *interp_view, *gamma_view;
	VecGetArray(interp, &interp_view);
	VecGetArray(gamma, &gamma_view);
	while (!blocks.empty()) {
		num_types++;
		// the first in the set is the type of interface that we are going to solve for
		set<Block> todo;
		Block      curr_type = *blocks.begin();
		blocks.erase(blocks.begin());
		todo.insert(curr_type);
		set<Block> to_be_deleted;
		for (auto iter = blocks.begin(); iter != blocks.end(); iter++) {
			if (*iter == curr_type) {
				todo.insert(*iter);
				to_be_deleted.insert(*iter);
			}
		}
		for (Block i : to_be_deleted) {
			blocks.erase(i);
		}

		// create domain representing curr_type
		Domain ds;
		ds.n         = n;
		ds.id        = 0;
		ds.id_local  = 0;
		ds.x_length  = 1;
		ds.y_length  = 1;
		ds.z_length  = 1;
		ds.nbr_id[0] = 1;
		ds.neumann   = curr_type.neumann;
		ds.local_i   = {{0, 0, 0, 0, 0, 0}};
		solver->addDomain(ds);

		map<BlockKey, shared_ptr<valarray<double>>> coeffs;
		// allocate blocks of coefficients
		for (const Block &b : todo) {
			shared_ptr<valarray<double>> ptr = coeffs[b];
			if (ptr.get() == nullptr) {
				coeffs[b] = shared_ptr<valarray<double>>(new valarray<double>(n * n * n * n));
			}
		}

		for (int j = 0; j < n * n; j++) {
			gamma_view[j] = 1;
			solver->solve(ds, f, u, gamma);
			gamma_view[j] = 0;

			// fill the blocks
			for (auto &p : coeffs) {
				Side      s    = p.first.s;
				IfaceType type = IfaceType::normal;
				VecScale(interp, 0);
				interpolator->interpolate(ds, s, type, u, interp);
				valarray<double> &block = *p.second;
				for (int i = 0; i < n * n; i++) {
					block[i * n * n + j] = -interp_view[i];
				}
				if (s == Side::west) {
					switch (type) {
						case IfaceType::normal:
							block[n * n * j + j] += 0.5;
							break;
						default:
							break;
					}
				}
			}
		}

		// now insert these results into the matrix for each interface
		for (Block block : todo) {
			insertBlock(&block, coeffs[block]);
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
	int     n = dc.n;
	PW<Mat> A;
	MatCreate(MPI_COMM_WORLD, &A);
	int local_size  = dc.ifaces.size() * n * n;
	int global_size = dc.num_global_interfaces * n * n;
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 19 * n * n, nullptr, 19 * n * n, nullptr);

	auto insertBlock = [&](Block *b, shared_ptr<valarray<double>> coeffs) {
		int global_i = b->i * n * n;
		int global_j = b->j * n * n;

		valarray<double> &orig = *coeffs;
		const function<int(int, int)> transforms_left[4]
		= {[&](int xi, int yi) { return xi + yi * n; },
		   [&](int xi, int yi) { return n - yi - 1 + xi * n; },
		   [&](int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
		   [&](int xi, int yi) { return yi + (n - xi - 1) * n; }};

		const function<int(int, int)> transforms_right[4]
		= {[&](int xi, int yi) { return xi + yi * n; },
		   [&](int xi, int yi) { return yi + (n - xi - 1) * n; },
		   [&](int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
		   [&](int xi, int yi) { return n - yi - 1 + xi * n; }};

		const function<int(int, int)> transforms_left_inv[4]
		= {[&](int xi, int yi) {
			   xi = n - xi - 1;
			   return xi + yi * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return n - yi - 1 + xi * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return n - xi - 1 + (n - yi - 1) * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return yi + (n - xi - 1) * n;
		   }};
		const function<int(int, int)> transforms_right_inv[4]
		= {[&](int xi, int yi) {
			   xi = n - xi - 1;
			   return xi + yi * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return yi + (n - xi - 1) * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return n - xi - 1 + (n - yi - 1) * n;
		   },
		   [&](int xi, int yi) {
			   xi = n - xi - 1;
			   return n - yi - 1 + xi * n;
		   }};

		valarray<double> copy(n * n * n * n);
		function<int(int, int)> col_trans, row_trans;
		if (b->mainLeft()) {
			if (b->mainFlipped()) {
				col_trans = transforms_left_inv[b->mainRot()];
			} else {
				col_trans = transforms_left[b->mainRot()];
			}
		} else {
			if (b->mainFlipped()) {
				col_trans = transforms_right_inv[b->mainRot()];
			} else {
				col_trans = transforms_right[b->mainRot()];
			}
		}
		if (b->auxLeft()) {
			if (b->auxFlipped()) {
				row_trans = transforms_left_inv[b->auxRot()];
			} else {
				row_trans = transforms_left[b->auxRot()];
			}
		} else {
			if (b->auxFlipped()) {
				row_trans = transforms_right_inv[b->auxRot()];
			} else {
				row_trans = transforms_right[b->auxRot()];
			}
		}
		for (int row_yi = 0; row_yi < n; row_yi++) {
			for (int row_xi = 0; row_xi < n; row_xi++) {
				int i_dest = row_xi + row_yi * n;
				int i_orig = row_trans(row_xi, row_yi);
				for (int col_yi = 0; col_yi < n; col_yi++) {
					for (int col_xi = 0; col_xi < n; col_xi++) {
						int j_dest = col_xi + col_yi * n;
						int j_orig = col_trans(col_xi, col_yi);

						copy[i_dest * n * n + j_dest] = orig[i_orig * n * n + j_orig];
					}
				}
			}
		}
		vector<int> inds_i(n * n);
		iota(inds_i.begin(), inds_i.end(), global_i);
		vector<int> inds_j(n * n);
		iota(inds_j.begin(), inds_j.end(), global_j);

		MatSetValues(A, n * n, &inds_i[0], n * n, &inds_j[0], &copy[0], ADD_VALUES);
	};

	assembleMatrix(insertBlock);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	return A;
}
PW_explicit<Mat> SchurHelper::formPBMatrix()
{
	int n           = dc.n;
	int local_size  = dc.ifaces.size() * n * n;
	int global_size = dc.num_global_interfaces * n * n;

	PBMatrix *APB    = new PBMatrix(n, local_size, global_size);
	auto insertBlock = [&](Block *b, shared_ptr<valarray<double>> coeffs) {
		int global_i = b->i;
		int global_j = b->j;

		const function<int(int, int, int)> transforms_left[4]
		= {[](int n, int xi, int yi) { return xi + yi * n; },
		   [](int n, int xi, int yi) { return n - yi - 1 + xi * n; },
		   [](int n, int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
		   [](int n, int xi, int yi) { return yi + (n - xi - 1) * n; }};

		const function<int(int, int, int)> transforms_right[4]
		= {[](int n, int xi, int yi) { return xi + yi * n; },
		   [](int n, int xi, int yi) { return yi + (n - xi - 1) * n; },
		   [](int n, int xi, int yi) { return n - xi - 1 + (n - yi - 1) * n; },
		   [](int n, int xi, int yi) { return n - yi - 1 + xi * n; }};

		const function<int(int, int, int)> transforms_left_inv[4]
		= {[](int n, int xi, int yi) {
			   xi = n - xi - 1;
			   return xi + yi * n;
		   },
		   [](int n, int xi, int yi) {
			   xi = n - xi - 1;
			   return n - yi - 1 + xi * n;
		   },
		   [](int n, int xi, int yi) {
			   xi = n - xi - 1;
			   return n - xi - 1 + (n - yi - 1) * n;
		   },
		   [](int n, int xi, int yi) {
			   xi = n - xi - 1;
			   return yi + (n - xi - 1) * n;
		   }};
		const function<int(int, int, int)> transforms_right_inv[4]
		= {[](int n, int xi, int yi) {
			   xi = n - xi - 1;
			   return xi + yi * n;
		   },
		   [](int n, int xi, int yi) {
			   xi = n - xi - 1;
			   return yi + (n - xi - 1) * n;
		   },
		   [](int n, int xi, int yi) {
			   xi = n - xi - 1;
			   return n - xi - 1 + (n - yi - 1) * n;
		   },
		   [](int n, int xi, int yi) {
			   xi = n - xi - 1;
			   return n - yi - 1 + xi * n;
		   }};

		function<int(int, int, int)> col_trans, row_trans;
		if (b->mainLeft()) {
			if (b->mainFlipped()) {
				col_trans = transforms_left_inv[b->mainRot()];
			} else {
				col_trans = transforms_left[b->mainRot()];
			}
		} else {
			if (b->mainFlipped()) {
				col_trans = transforms_right_inv[b->mainRot()];
			} else {
				col_trans = transforms_right[b->mainRot()];
			}
		}
		if (b->auxLeft()) {
			if (b->auxFlipped()) {
				row_trans = transforms_left_inv[b->auxRot()];
			} else {
				row_trans = transforms_left[b->auxRot()];
			}
		} else {
			if (b->auxFlipped()) {
				row_trans = transforms_right_inv[b->auxRot()];
			} else {
				row_trans = transforms_right[b->auxRot()];
			}
		}
		APB->insertBlock(global_i, global_j, coeffs, col_trans, row_trans);
	};

	assembleMatrix(insertBlock);
	APB->finalize();

	return APB->getMatrix();
}
