#include "SchurHelper.h"
#include <array>
#include <iostream>
#include <numeric>
#include <petscao.h>
#include <tuple>
using namespace std;

SchurHelper::SchurHelper(DomainCollection dc, shared_ptr<PatchSolver> solver,
                         shared_ptr<PatchOperator> op, shared_ptr<Interpolator> interpolator)
{
	this->n = dc.getN();
	for (auto &p : dc.domains) {
		domains.push_back(*p.second);
	}
	map<int, pair<int, IfaceSet>> off_proc_ifaces;
	for (SchurDomain &sd : domains) {
		sd.enumerateIfaces(ifaces, off_proc_ifaces);
		solver->addDomain(sd);
	}
	{
		// send info
		deque<char *>       buffers;
		deque<char *>       recv_buffers;
		vector<MPI_Request> requests;
		for (auto &p : off_proc_ifaces) {
			int       dest   = p.second.first;
			IfaceSet &iface  = p.second.second;
			int       size   = iface.serialize(nullptr);
			char *    buffer = new char[size];
			buffers.push_back(buffer);
			iface.serialize(buffer);
			MPI_Request request;
			MPI_Isend(buffer, size, MPI_CHAR, dest, 0, MPI_COMM_WORLD, &request);
			requests.push_back(request);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		int        is_message;
		MPI_Status status;
		MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &is_message, &status);
		// recv info
		while (is_message) {
			int size;
			MPI_Get_count(&status, MPI_CHAR, &size);
			char *buffer = new char[size];
			recv_buffers.push_back(buffer);

			MPI_Request request;
			MPI_Irecv(buffer, size, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
			          &request);
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &is_message, &status);
			requests.push_back(request);
		}
		// wait for all
		vector<MPI_Status> statuses(requests.size());
		MPI_Waitall(requests.size(), &requests[0], &statuses[0]);
		// delete send buffers
		for (char *buffer : buffers) {
			delete[] buffer;
		}
		// process received objects
		for (char *buffer : recv_buffers) {
			IfaceSet ifs = IfaceSet::deserialize(buffer);
			ifaces[ifs.id].insert(ifs);
			delete[] buffer;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	indexDomainIfacesLocal();
	indexIfacesLocal();
	this->solver       = solver;
	this->op           = op;
	this->interpolator = interpolator;
	local_gamma        = getNewSchurDistVec();
	local_interp       = getNewSchurDistVec();
	gamma              = getNewSchurVec();
	PW<IS> dist_is;
	ISCreateBlock(MPI_COMM_SELF, n * n, iface_dist_map_vec.size(), &iface_dist_map_vec[0],
	              PETSC_COPY_VALUES, &dist_is);
	VecScatterCreate(gamma, dist_is, local_gamma, nullptr, &scatter);

	int num_ifaces = ifaces.size();
	MPI_Allreduce(&num_ifaces, &num_global_ifaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	cerr << "Num Global All: " << num_global_ifaces << endl;
	cerr << "Num Local All: " << num_ifaces << endl;
}

void SchurHelper::zoltanBalance() {}
void SchurHelper::solveWithInterface(const Vec f, Vec u, const Vec gamma, Vec diff)
{
	// initilize our local variables
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	VecScale(local_interp, 0);
	// solve over domains on this proc
	for (SchurDomain &sd : domains) {
		solver->solve(sd, f, u, local_gamma);
		interpolator->interpolate(sd, u, local_interp);
	}

	// export diff vector
	VecScale(diff, 0);
	VecScatterBegin(scatter, local_interp, diff, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, diff, ADD_VALUES, SCATTER_REVERSE);
	VecAXPBY(diff, 1.0, -1.0, gamma);
}
void SchurHelper::solveAndInterpolateWithInterface(const Vec f, Vec u, const Vec gamma, Vec interp)
{
	// initilize our local variables
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	VecScale(local_interp, 0);
	// solve over domains on this proc
	for (SchurDomain &sd : domains) {
		solver->solve(sd, f, u, local_gamma);
		interpolator->interpolate(sd, u, local_interp);
	}

	// export diff vector
	VecScale(interp, 0);
	VecScatterBegin(scatter, local_interp, interp, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, interp, ADD_VALUES, SCATTER_REVERSE);
}
void SchurHelper::solveWithSolution(const Vec f, Vec u)
{
	// initilize our local variables
	VecScale(local_gamma, 0);
	VecScale(local_interp, 0);
	for (SchurDomain &sd : domains) {
		interpolator->interpolate(sd, u, local_interp);
	}
	VecScale(gamma, 0);
	VecScatterBegin(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	// solve over domains on this proc
	for (SchurDomain &sd : domains) {
		solver->solve(sd, f, u, local_gamma);
	}
}
void SchurHelper::interpolateToInterface(const Vec f, Vec u, Vec gamma)
{
	// initilize our local variables
	VecScale(local_interp, 0);
	for (SchurDomain &sd : domains) {
		interpolator->interpolate(sd, u, local_interp);
	}
	VecScale(gamma, 0);
	VecScatterBegin(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
}
void SchurHelper::applyWithInterface(const Vec u, const Vec gamma, Vec f)
{
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	for (SchurDomain &sd : domains) {
		op->apply(sd, u, local_gamma, f);
	}
}
void SchurHelper::apply(const Vec u, Vec f)
{
	VecScale(local_interp, 0);
	for (SchurDomain &sd : domains) {
		interpolator->interpolate(sd, u, local_interp);
	}
	VecScale(gamma, 0);
	VecScatterBegin(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterEnd(scatter, local_interp, gamma, ADD_VALUES, SCATTER_REVERSE);
	VecScatterBegin(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, gamma, local_gamma, INSERT_VALUES, SCATTER_FORWARD);

	for (SchurDomain &sd : domains) {
		op->apply(sd, u, local_gamma, f);
	}
}
enum class Rotation : char { x_cw, x_ccw, y_cw, y_ccw, z_cw, z_ccw };
struct Block {
	static const Side             side_table[6][6];
	static const char             rots_table[6][6];
	static const vector<Rotation> main_rot_plan[6];
	static const vector<Rotation> aux_rot_plan_dirichlet[6];
	static const vector<Rotation> aux_rot_plan_neumann[16];
	static const char             rot_quad_lookup_left[4][4];
	static const char             rot_quad_lookup_right[4][4];
	static const char             quad_flip_lookup[4];
	char                          data = 0;
	IfaceType                     type;
	Side                          main;
	Side                          aux;
	int                           j;
	int                           i;
	bitset<6>                     neumann;
	Block(Side main, int j, Side aux, int i, bitset<6> neumann, IfaceType type)
	{
		this->main    = main;
		this->j       = j;
		this->aux     = aux;
		this->i       = i;
		this->neumann = neumann;
		this->type    = type;
		data          = 0;
		data |= main.isLowerOnAxis() << 7;
		data |= aux.isLowerOnAxis() << 3;
		rotate();
	}
	void applyRotation(const Rotation rot)
	{
		char r;
		// main rotation
		r = (data & ~(~0u << 2) << 4) >> 4;
		r = (r + rots_table[static_cast<int>(rot)][main.toInt()]) & 0b11;
		data &= ~(~(~0u << 2) << 4);
		data |= r << 4;
		// aux rotation
		r = data & ~(~0u << 2);
		r = (r + rots_table[static_cast<int>(rot)][aux.toInt()]) & 0b11;
		data &= ~0u << 2;
		data |= r;
		main                  = side_table[static_cast<int>(rot)][main.toInt()];
		aux                   = side_table[static_cast<int>(rot)][aux.toInt()];
		bitset<6> old_neumann = neumann;
		for (int i = 0; i < 6; i++) {
			neumann[side_table[(int) rot][i].toInt()] = old_neumann[i];
		}
	}
	void rotate()
	{
		for (Rotation rot : main_rot_plan[main.toInt()]) {
			applyRotation(rot);
		}
		if (neumann.to_ulong() == 0) {
			for (Rotation rot : aux_rot_plan_dirichlet[aux.toInt()]) {
				applyRotation(rot);
			}
		} else {
			for (Rotation rot : aux_rot_plan_neumann[neumann.to_ulong() >> 2]) {
				applyRotation(rot);
			}
		}
		// updated iface type
		auto rotateQuad = [&](int quad) {
			if (auxOrigLeft()) {
				quad = rot_quad_lookup_left[auxRot()][quad];
			} else {
				quad = rot_quad_lookup_right[auxRot()][quad];
			}
			if (auxFlipped()) { quad = quad_flip_lookup[quad]; }
			return quad;
		};
		switch (type) {
			case IfaceType::fine_to_coarse_0:
			case IfaceType::fine_to_coarse_1:
			case IfaceType::fine_to_coarse_2:
			case IfaceType::fine_to_coarse_3: {
				int quad = (int) type - (int) IfaceType::fine_to_coarse_0;
				quad     = rotateQuad(quad);
				type     = IfaceType::fine_to_coarse_0 + quad;
			} break;
			case IfaceType::fine_to_fine_0:
			case IfaceType::fine_to_fine_1:
			case IfaceType::fine_to_fine_2:
			case IfaceType::fine_to_fine_3: {
				int quad = (int) type - (int) IfaceType::fine_to_fine_0;
				quad     = rotateQuad(quad);
				type     = IfaceType::fine_to_fine_0 + quad;
			} break;
			case IfaceType::coarse_to_fine_0:
			case IfaceType::coarse_to_fine_1:
			case IfaceType::coarse_to_fine_2:
			case IfaceType::coarse_to_fine_3: {
				int quad = (int) type - (int) IfaceType::coarse_to_fine_0;
				quad     = rotateQuad(quad);
				type     = IfaceType::coarse_to_fine_0 + quad;
			} break;
			default:
				break;
		}
	}
	bool operator<(const Block &b) const
	{
		return std::tie(i, j, data) < std::tie(b.i, b.j, b.data);
	}
	bool operator==(const Block &b) const
	{
		return neumann.to_ulong() == b.neumann.to_ulong();
	}
	//
	bool mainLeft()
	{
		return main.isLowerOnAxis();
	}
	bool mainFlipped()
	{
		bool left      = main.isLowerOnAxis();
		bool orig_left = (data >> 7) & 0b1;
		return left != orig_left;
	}
	int mainRot()
	{
		return (data >> 4) & 0b11;
	}
	bool auxOrigLeft()
	{
		return (data >> 3) & 0b1;
	}
	bool auxLeft()
	{
		return aux.isLowerOnAxis();
	}
	bool auxFlipped()
	{
		bool left      = aux.isLowerOnAxis();
		bool orig_left = (data >> 3) & 0b1;
		return left != orig_left;
	}
	int auxRot()
	{
		return data & 0b11;
	}
};
struct BlockKey {
	IfaceType type;
	Side      s;

	BlockKey() {}
	BlockKey(const Block &b)
	{
		type = b.type;
		s    = b.aux;
	}
	friend bool operator<(const BlockKey &l, const BlockKey &r)
	{
		return tie(l.s, l.type) < tie(r.s, r.type);
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
const vector<Rotation> Block::aux_rot_plan_dirichlet[6]
= {{}, {}, {}, {Rotation::x_cw, Rotation::x_cw}, {Rotation::x_cw}, {Rotation::x_ccw}};
const vector<Rotation> Block::aux_rot_plan_neumann[16] = {{},
                                                          {},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {Rotation::x_cw, Rotation::x_cw},
                                                          {Rotation::x_ccw},
                                                          {},
                                                          {Rotation::x_ccw},
                                                          {},
                                                          {Rotation::x_ccw},
                                                          {Rotation::x_cw},
                                                          {Rotation::x_ccw},
                                                          {}};
const char             Block::rot_quad_lookup_left[4][4]
= {{0, 1, 2, 3}, {1, 3, 0, 2}, {3, 2, 1, 0}, {2, 0, 3, 1}};
const char Block::rot_quad_lookup_right[4][4]
= {{0, 1, 2, 3}, {2, 0, 3, 1}, {3, 2, 1, 0}, {1, 3, 0, 2}};
const char Block::quad_flip_lookup[4] = {1, 0, 3, 2};
void       SchurHelper::assembleMatrix(inserter insertBlock)
{
	set<Block> blocks;
	for (auto &p : ifaces) {
		IfaceSet &ifs = p.second;
		int       i   = ifs.id_global;
		for (const Iface &iface : ifs.ifaces) {
			Side aux = iface.s;
			for (int s = 0; s < 6; s++) {
				int j = iface.global_id[s];
				if (j != -1) {
					Side main = static_cast<Side>(s);
					blocks.insert(Block(main, j, aux, i, iface.neumann, iface.type));
				}
			}
		}
	}
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
		SchurDomain sd;
		sd.n                           = n;
		sd.x_length                    = 1;
		sd.y_length                    = 1;
		sd.z_length                    = 1;
		sd.neumann                     = curr_type.neumann;
		sd.getIfaceInfoPtr(Side::west) = new NormalIfaceInfo();
		solver->addDomain(sd);

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
			solver->solve(sd, f, u, gamma);
			gamma_view[j] = 0;

			// fill the blocks
			for (auto &p : coeffs) {
				BlockKey  bk   = p.first;
				Side      s    = bk.s;
				IfaceType type = bk.type;
				VecScale(interp, 0);
				interpolator->interpolate(sd, s, 0, type, u, interp);
				valarray<double> &block = *p.second;
				for (int i = 0; i < n * n; i++) {
					block[i * n * n + j] = -interp_view[i];
				}
				if (s == Side::west) {
					switch (type) {
						case IfaceType::normal:
							block[n * n * j + j] += 0.5;
							break;
						case IfaceType::coarse_to_coarse:
						case IfaceType::fine_to_fine_0:
						case IfaceType::fine_to_fine_1:
						case IfaceType::fine_to_fine_2:
						case IfaceType::fine_to_fine_3:
							block[n * n * j + j] += 1.0;
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
	PW<Mat> A;
	MatCreate(MPI_COMM_WORLD, &A);
	int local_size  = ifaces.size() * n * n;
	int global_size = getSchurVecGlobalSize();
	MatSetSizes(A, local_size, local_size, global_size, global_size);
	MatSetType(A, MATMPIAIJ);
	MatMPIAIJSetPreallocation(A, 10 * n * n, nullptr, 10 * n * n, nullptr);

	auto insertBlock = [&](Block *b, shared_ptr<valarray<double>> coeffs) {
		int global_i = b->i * n * n;
		int global_j = b->j * n * n;

		valarray<double> &            orig = *coeffs;
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

		valarray<double>        copy(n * n * n * n);
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
PBMatrix *SchurHelper::formPBMatrix()
{
	int local_size  = ifaces.size() * n * n;
	int global_size = getSchurVecGlobalSize();

	PBMatrix *APB         = new PBMatrix(n, local_size, global_size);
	auto      insertBlock = [&](Block *b, shared_ptr<valarray<double>> coeffs) {
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

	return APB;
}
PW_explicit<Mat> SchurHelper::getPBMatrix()
{
	return formPBMatrix()->getMatrix();
}
PW_explicit<Mat> SchurHelper::getPBDiagInv()
{
	return formPBMatrix()->getDiagInv()->getMatrix();
}
void SchurHelper::getPBDiagInv(PC p)
{
	formPBMatrix()->getDiagInv()->getPrec(p);
}
void SchurHelper::indexDomainIfacesLocal()
{
	vector<int>   map_vec;
	map<int, int> rev_map;
	if (!domains.empty()) {
		int curr_i = 0;
		for (SchurDomain &sd : domains) {
			for (int id : sd.getIds()) {
				if (rev_map.count(id) == 0) {
					rev_map[id] = curr_i;
					map_vec.push_back(id);
					curr_i++;
				}
			}
		}
		for (SchurDomain &sd : domains) {
			sd.setLocalIndexes(rev_map);
		}
	}
	iface_dist_map_vec = map_vec;
}
void SchurHelper::indexIfacesLocal()
{
	int           curr_i = 0;
	vector<int>   map_vec;
	vector<int>   off_proc_map_vec;
	vector<int>   off_proc_map_vec_send;
	map<int, int> rev_map;
	if (!ifaces.empty()) {
		set<int> todo;
		for (auto &p : ifaces) {
			todo.insert(p.first);
		}
		set<int> enqueued;
		while (!todo.empty()) {
			deque<int> queue;
			queue.push_back(*todo.begin());
			enqueued.insert(*todo.begin());
			while (!queue.empty()) {
				int i = queue.front();
				todo.erase(i);
				queue.pop_front();
				map_vec.push_back(i);
				IfaceSet &ifs = ifaces.at(i);
				rev_map[i]    = curr_i;
				curr_i++;
				for (int nbr : ifs.getNbrs()) {
					if (!enqueued.count(nbr)) {
						enqueued.insert(nbr);
						if (ifaces.count(nbr)) {
							queue.push_back(nbr);
						} else {
							off_proc_map_vec.push_back(nbr);
						}
					}
				}
			}
		}
	}

	// map off proc
	for (int i : off_proc_map_vec) {
		rev_map[i] = curr_i;
		curr_i++;
	}
	for (auto &p : ifaces) {
		p.second.setLocalIndexes(rev_map);
	}
	iface_map_vec          = map_vec;
	iface_off_proc_map_vec = off_proc_map_vec;
	indexIfacesGlobal();
}
void SchurHelper::indexIfacesGlobal()
{
	// global indices are going to be sequentially increasing with rank
	int local_size = ifaces.size();
	int start_i;
	MPI_Scan(&local_size, &start_i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	start_i -= local_size;
	vector<int> new_global(local_size);
	iota(new_global.begin(), new_global.end(), start_i);

	// create map for gids
	PW<AO> ao;
	AOCreateMapping(MPI_COMM_WORLD, local_size, &iface_map_vec[0], &new_global[0], &ao);

	// get indices for schur matrix
	{
		// get global indices that we want to recieve for dest vector
		vector<int> inds = iface_map_vec;
		for (int i : iface_off_proc_map_vec) {
			inds.push_back(i);
		}

		// get new global indices
		AOApplicationToPetsc(ao, inds.size(), &inds[0]);
		map<int, int> rev_map;
		for (size_t i = 0; i < inds.size(); i++) {
			rev_map[i] = inds[i];
		}

		// set new global indices in iface objects
		for (auto &p : ifaces) {
			p.second.setGlobalIndexes(rev_map);
		}
		for (size_t i = 0; i < iface_map_vec.size(); i++) {
			iface_map_vec[i] = inds[i];
		}
		for (size_t i = 0; i < iface_off_proc_map_vec.size(); i++) {
			iface_off_proc_map_vec[i] = inds[iface_map_vec.size() + i];
		}
	}
	// get indices for local ifaces
	{
		// get global indices that we want to recieve for dest vector
		vector<int> inds = iface_dist_map_vec;

		// get new global indices
		AOApplicationToPetsc(ao, inds.size(), &inds[0]);
		map<int, int> rev_map;
		for (size_t i = 0; i < inds.size(); i++) {
			rev_map[i] = inds[i];
		}

		// set new global indices in domain objects
		for (SchurDomain &sd : domains) {
			sd.setGlobalIndexes(rev_map);
		}
		for (size_t i = 0; i < iface_dist_map_vec.size(); i++) {
			iface_dist_map_vec[i] = inds[i];
		}
	}
}
PW_explicit<Vec> SchurHelper::getNewSchurVec()
{
	PW<Vec> u;
	VecCreateMPI(MPI_COMM_WORLD, iface_map_vec.size() * n * n, PETSC_DETERMINE, &u);
	return u;
}
PW_explicit<Vec> SchurHelper::getNewSchurDistVec()
{
	PW<Vec> u;
	VecCreateSeq(PETSC_COMM_SELF, iface_dist_map_vec.size() * n * n, &u);
	return u;
}
