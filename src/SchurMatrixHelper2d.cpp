#include "SchurMatrixHelper2d.h"
#include <array>
#include <numeric>
#include <tuple>
using namespace std;
enum axis_enum { X_AXIS, Y_AXIS };
enum bc_enum { DIRICHLET, NEUMANN, REFINED };

struct Block {
	IfaceType           type;
	Side<2>             s;
	bitset<4>           neumann;
	int                 i;
	int                 j;
	bool                flip_i;
	bool                flip_j;
	bitset<4>           flip_j_table = 0b0110;
	array<bitset<4>, 4> flip_i_table = {{0b0000, 0b1111, 0b0011, 0b1100}};
	Side<2> rot_table[4][4] = {{Side<2>::west, Side<2>::east, Side<2>::south, Side<2>::north},
	                           {Side<2>::east, Side<2>::west, Side<2>::north, Side<2>::south},
	                           {Side<2>::north, Side<2>::south, Side<2>::west, Side<2>::east},
	                           {Side<2>::south, Side<2>::north, Side<2>::east, Side<2>::west}};

	Block(Side<2> main, int j, Side<2> aux, int i, bitset<4> neumann, IfaceType type)
	{
		s       = rot_table[main.toInt()][aux.toInt()];
		this->i = i;
		this->j = j;
        this->type = type;
        for(int s=0;s<4;s++){
            this->neumann[rot_table[main.toInt()][s].toInt()]=neumann[s];
        }
		flip_j  = flip_j_table[main.toInt()];
		flip_i  = flip_i_table[main.toInt()][s.toInt()];
        if(flip_i){
            this->type.setOrthant(!type.getOrthant());
        }
	}
	bool operator==(const Block &b) const
	{
		return neumann.to_ulong() == b.neumann.to_ulong();
	}
	bool operator<(const Block &b) const
	{
		return std::tie(i, j, flip_j) < std::tie(b.i, b.j,b.flip_j);
	}
};
struct BlockKey {
	IfaceType type;
	Side<2>   s;

	BlockKey() {}
	BlockKey(const Block &b)
	{
		s    = b.s;
		type = b.type;
	}
	friend bool operator<(const BlockKey &l, const BlockKey &r)
	{
		return std::tie(l.s, l.type) < std::tie(r.s, r.type);
	}
};
void SchurMatrixHelper2d::assembleMatrix(inserter insertBlock)
{
	set<Block> blocks;
	for (auto p : sh->getIfaces()) {
		const IfaceSet<2> &ifs = p.second;
		int                i   = ifs.id_global;
		for (const Iface<2> &iface : ifs.ifaces) {
			Side<2> aux = iface.s;
			for (int s = 0; s < 4; s++) {
				int j = iface.global_id[s];
				if (j != -1) {
					Side<2> main = static_cast<Side<2>>(s);
					blocks.insert(Block(main, j, aux, i, iface.neumann, iface.type));
				}
			}
		}
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

		auto solver       = sh->getSolver();
		auto interpolator = sh->getInterpolator();
		// create domain representing curr_type
		SchurDomain<2> sd;
		sd.n = n;
		sd.domain.n = n;
		sd.neumann                        = curr_type.neumann;
		sd.getIfaceInfoPtr(Side<2>::west) = new NormalIfaceInfo<2>();
		solver->addDomain(sd);
		std::deque<SchurDomain<2>> single_domain;
		single_domain.push_back(sd);

		map<BlockKey, shared_ptr<valarray<double>>> coeffs;
		// allocate blocks of coefficients
		for (const Block &b : todo) {
			shared_ptr<valarray<double>> ptr = coeffs[b];
			if (ptr.get() == nullptr) {
				coeffs[b] = shared_ptr<valarray<double>>(new valarray<double>(n * n));
			}
		}

		for (int j = 0; j < n; j++) {
			gamma_view[j] = 1;
			solver->domainSolve(single_domain, f, u, gamma);
			gamma_view[j] = 0;

			// fill the blocks
			for (auto &p : coeffs) {
				Side<2>   s    = p.first.s;
				IfaceType type = p.first.type;
				VecScale(interp, 0);
				interpolator->interpolate(sd, s, 0, type, u, interp);
				valarray<double> &block = *p.second;
				for (int i = 0; i < n; i++) {
					block[i * n + j] = -interp_view[i];
				}
				if (s == Side<2>::west) {
					switch (type.toInt()) {
						case IfaceType::normal:
							block[n * j + j] += 0.5;
							break;
						case IfaceType::coarse_to_coarse:
						case IfaceType::fine_to_fine:
							block[n * j + j] += 1;
							break;
						default:
							break;
					}
				}
			}
		}

		// now insert these results into the matrix for each interface
		for (Block block : todo) {
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
PW_explicit<Mat> SchurMatrixHelper2d::formCRSMatrix()
{
	PW<Mat> A;
	MatCreate(MPI_COMM_WORLD, &A);
	int local_size  = sh->getIfaces().size() * n;
	int global_size = sh->getSchurVecGlobalSize();
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
			  if (flip_i) { block_i = n - i - 1; }
			  for (int j = 0; j < n; j++) {
				  int block_j = j;
				  if (flip_j) { block_j = n - j - 1; }
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
