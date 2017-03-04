#include "RBMatrix.h"


extern "C" {

// LU decomoposition of a general matrix
void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *lda, int *IPIV, double *B, int *LDB,
             int *INFO);
void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA,
            double *B, int *LDB, double *BETA, double *C, int *LDC);

// generate inverse of a matrix given its LU decomposition
void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork, int *INFO);
}
using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
typedef RCP<valarray<double>> blk_ptr;
RBMatrix::RBMatrix(Teuchos::RCP<map_type> map, int block_size, int num_blocks)
{
	domain = map;
	range  = map;

	this->num_blocks = num_blocks;
	this->block_size = block_size;

	block_cols = vector<std::map<int, Block>>(num_blocks);
	for (int x = 0; x < num_blocks; x++) {
		int i = domain->getGlobalElement(curr_local_i);
		range_map[i] = curr_local_i;
		for (int x = i; x < i + block_size; x++) {
			global_i.push_back(x);
		}
		curr_local_i += block_size;
	}
}
void LUSolver::apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode, double alpha,
                     double beta) const
{
	auto x_view = x.getLocalView<Kokkos::HostSpace>();
	auto y_view = y.getLocalView<Kokkos::HostSpace>();
	int  block_size = L->getBlockSize();
	//L solve
	for (int j = 0; j < L->getNumBlocks(); j++) {
		// do lower triangular block
		blk_ptr Lkk = L->getBlock(j * block_size, j * block_size);
        int x_i = x.getMap()->getLocalElement(j*block_size);

		for (int block_i = 0; block_i < block_size; block_i++) {
			for (int block_j = block_i - 1; block_j >= 0; block_j--) {
				x_view(x_i + block_i, 0)
				-= x_view(x_i + block_j, 0) * (*Lkk)[block_i + block_size * block_j];
			}
		}

		// do rest of blocks
		for (int i = j + 1; i < L->getNumBlocks(); i++) {
			int     x_i2 = x.getMap()->getLocalElement(i * block_size);
			blk_ptr Lik  = L->getBlock(i * block_size, j * block_size);
			if (!Lik.is_null()) {
				for (int block_j = 0; block_j < block_size; block_j++) {
					for (int block_i = 0; block_i < block_size; block_i++) {
						x_view(x_i2 + block_i, 0)
						-= x_view(x_i + block_j, 0) * (*Lik)[block_i + block_size * block_j];
					}
				}
			}
		}
	}
	// U solve
	for (int j = U->getNumBlocks() - 1; j >= 0; j--) {
		// do lower triangular block
		blk_ptr Ukk = U->getBlock(j * block_size, j * block_size);
		int     x_i = x.getMap()->getLocalElement(j * block_size);

		for (int block_i = block_size - 1; block_i >= 0; block_i--) {
			for (int block_j = block_size - 1; block_j > block_i; block_j--) {
				x_view(x_i + block_i, 0)
				-= x_view(x_i + block_j, 0) * (*Ukk)[block_i + block_size * block_j];
			}
			// dvide by diagonal
			x_view(x_i + block_i, 0) /= (*Ukk)[block_i + block_size * block_i];
		}

		// do rest of blocks
		for (int i = j - 1; i >= 0; i--) {
			int     x_i2 = x.getMap()->getLocalElement(i * block_size);
			blk_ptr Uik  = U->getBlock(i * block_size, j * block_size);
			if (!Uik.is_null()) {
				for (int block_j = 0; block_j < block_size; block_j++) {
					for (int block_i = 0; block_i < block_size; block_i++) {
						x_view(x_i2 + block_i, 0)
						-= x_view(x_i + block_j, 0) * (*Uik)[block_i + block_size * block_j];
						//cerr << x_view(x_i2 + block_i, 0) << endl;
					}
				}
			}
		}
	}
	//cerr << "hello\n";

	y.update(1,x,0);
}
void RBMatrix::apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode, double alpha,
                     double beta) const
{
    vector_type my_y(range, 1);
	auto x_view = x.getLocalView<Kokkos::HostSpace>();
	y.putScalar(0);
	my_y.putScalar(0);
	auto y_view = my_y.getLocalView<Kokkos::HostSpace>();
	// column loop
	for (int index = 0; index < num_blocks; index++) {
		int              start_j = index * block_size;
		valarray<double> slot_rev(block_size);
		// go down the column
		for (auto &p : block_cols[index]) {
			int   start_i    = p.first;
			Block curr_block = p.second;
			for (int iii = 0; iii < block_size; iii++) {
				int matrix_i = start_i + iii;
				if (matrix_i == local_skip_i && local_skip_i != -1 && local_skip_j != -1) {
					y_view(matrix_i, 0) += x_view(local_skip_j, 0);
				} else {
					int block_i = iii;
					if (curr_block.flip_i) {
						block_i = block_size - 1 - iii;
					}
					valarray<double> &curr_blk = *curr_block.block;

					if (curr_block.flip_j) {
						for (int i = 0; i < block_size; i++) {
							y_view(matrix_i, 0)
							+= x_view(start_j + i, 0)
							   * curr_blk[block_i + block_size * (block_size - 1 - i)];
						}
					} else {
						for (int i = 0; i < block_size; i++) {
							y_view(matrix_i, 0)
							+= x_view(start_j + i, 0) * curr_blk[block_i + block_size * i];
						}
					}
				}
			}
		}
	}
	y.doExport(my_y, *exporter, Tpetra::CombineMode::ADD);
}
void RBMatrix::insertBlock(int i, int j, RCP<valarray<double>> block, bool flip_i, bool flip_j)
{
	int   local_j = domain->getLocalElement(j);
	int   local_i = -1;
	try{
		local_i = range_map.at(i);
	} catch (const out_of_range &oor) {
		// do nothing
		// cerr << i << " is mapped to " << curr_local_i << "\n";
		local_i = curr_local_i;
		for (int x = i; x <i+ block_size; x++) {
			global_i.push_back(x);
		}
        curr_local_i+=block_size;
		range_map[i] = local_i;
	}

	Block b(block, flip_i, flip_j);

    int index = local_j/block_size;
	if (block_cols[index].count(local_i) > 0) {
		Block first_block = block_cols[index][local_i];

		pair<Block, Block> bpair(first_block, b);
		pair<Block, Block> bpair_rev(b, first_block);

		if (combos.count(bpair) > 0) {
			block_cols[index][local_i] = combos[bpair];
		} else if (combos.count(bpair_rev) > 0) {
			block_cols[index][local_i] = combos[bpair_rev];
		} else {
			RCP<valarray<double>> blk_ptr = rcp(new valarray<double>(first_block.block->size()));
			Block                 new_block(blk_ptr, false, false);
			valarray<double> &    new_blk    = *blk_ptr;
			valarray<double> &    first_blk  = *first_block.block;
			valarray<double> &    second_blk = *b.block;

			for (int i = 0; i < block_size; i++) {
				int block_i = i;
				if (first_block.flip_i) {
					block_i = block_size - i - 1;
				}
				for (int j = 0; j < block_size; j++) {
					int block_j = j;
					if (first_block.flip_j) {
						block_j = block_size - j - 1;
					}
					new_blk[i + block_size * j] = first_blk[block_i + block_size * block_j];
				}
			}
			for (int i = 0; i < block_size; i++) {
				int block_i = i;
				if (b.flip_i) {
					block_i = block_size - i - 1;
				}
				for (int j = 0; j < block_size; j++) {
					int block_j = j;
					if (b.flip_j) {
						block_j = block_size - j - 1;
					}
					new_blk[i + block_size * j] += second_blk[block_i + block_size * block_j];
				}
			}
			block_cols[index][local_i] = new_block;
			combos[bpair]        = new_block;
		}
	} else {
        nz+=block->size();
		block_cols[index][local_i] = b;
	}
}

void RBMatrix::createRangeMap()
{
	range = Teuchos::rcp(new map_type(-1, &global_i[0], global_i.size(), 0, domain->getComm()));
	if (skip_index != -1) {
		local_skip_i = range->getLocalElement(skip_index);
		local_skip_j = domain->getLocalElement(skip_index);
	}
	exporter = rcp(new Tpetra::Export<>(range, domain));
}
RCP<RBMatrix> RBMatrix::invBlockDiag(){
	RCP<RBMatrix>    Inv    = rcp(new RBMatrix(domain, block_size, num_blocks));
	map<double *, RCP<valarray<double>>> computed;
	for (int index = 0; index < num_blocks; index++) {
		int              start_j  = index * block_size;
		int              global_j = domain->getGlobalElement(start_j);
		// go down the column
		for (auto &p : block_cols[index]) {
			int   start_i    = p.first;
			int   global_i   = range->getGlobalElement(start_i);
			if (global_j == global_i) {
				RCP<valarray<double>> blk_inv_ptr;
				Block                 curr_block = p.second;
				if (local_skip_i != -1 && local_skip_i <= start_i
				    && start_i < local_skip_i + block_size) {
					// modify the row of this block and invert
					blk_inv_ptr               = rcp(new valarray<double>(*curr_block.block));
					valarray<double> &blk_inv = *blk_inv_ptr;

					// modify the row
					int block_i = local_skip_i % block_size;
					for (int j = 0; j < block_size; j++) {
						blk_inv[block_i * block_size + j] = 0;
					}
					blk_inv[block_i * block_size + block_i] = 1;

					// compute inverse of block
					valarray<int>    ipiv(block_size + 1);
					int              lwork = block_size * block_size;
					valarray<double> work(lwork);
					int              info;
					dgetrf_(&block_size, &block_size, &blk_inv[0], &block_size, &ipiv[0], &info);
					dgetri_(&block_size, &blk_inv[0], &block_size, &ipiv[0], &work[0], &lwork,
					        &info);
				} else {
					try {
						blk_inv_ptr = computed.at(&(*curr_block.block)[0]);
					} catch (const out_of_range &oor) {
						// invert this block
						blk_inv_ptr               = rcp(new valarray<double>(*curr_block.block));
						valarray<double> &blk_inv = *blk_inv_ptr;

						// compute inverse of block
						valarray<int>    ipiv(block_size + 1);
						int              lwork = block_size * block_size;
						valarray<double> work(lwork);
						int              info;
						dgetrf_(&block_size, &block_size, &blk_inv[0], &block_size, &ipiv[0],
						        &info);
						dgetri_(&block_size, &blk_inv[0], &block_size, &ipiv[0], &work[0], &lwork,
						        &info);

						// insert the block
						computed[&(*curr_block.block)[0]] = blk_inv_ptr;
					}
				}
				Inv->insertBlock(global_i, global_j, blk_inv_ptr, false, false);
			}
		}
	}
	Inv->createRangeMap();
	return Inv;
}
typedef RCP<valarray<int>> piv_ptr;

void RBMatrix::DGETRF(blk_ptr &A, blk_ptr &L, blk_ptr &U, piv_ptr &P)
{
	// compute inverse of block
	P         = rcp(new valarray<int>(block_size + 1));
	int lwork = block_size * block_size;

	valarray<double> work(lwork);
	int              info;
	L = rcp(new valarray<double>(*A));
	dgetrf_(&block_size, &block_size, &(*L)[0], &block_size, &(*P)[0], &info);
	U = rcp(new valarray<double>(*L));
	for (int i = 0; i < block_size; i++) {
		for (int j = i+1; j < block_size; j++) {
			(*L)[i + block_size * j] = 0;
		}
		(*L)[i + block_size * i] = 1;
	}
	for (int i = 0; i < block_size; i++) {
		for (int j = 0; j < i; j++) {
			(*U)[i + block_size * j] = 0;
		}
	}
}

void RBMatrix::DGESSM(blk_ptr &A, blk_ptr &L, blk_ptr &U, piv_ptr &P)
{
	U = rcp(new valarray<double>(*A));
	// L^-1*P*A
	char trans = 'N';
	int  info;
	dgetrs_(&trans, &block_size, &block_size, &(*L)[0], &block_size, &(*P)[0], &(*U)[0],
	        &block_size, &info);
}
void RBMatrix::LSolve(blk_ptr &A, blk_ptr &L, blk_ptr &U, piv_ptr &P)
{
	L = rcp(new valarray<double>(*A));
	valarray<double> TMP(*A);
	for (int j = 0; j < block_size; j++) {
		for (int i = 0; i < block_size; i++) {
			TMP[j + block_size * i] = (*A)[i + block_size * j];
		}
	}
	// L^-1*P*A
	char trans = 'T';
	int  info;
	dgetrs_(&trans, &block_size, &block_size, &(*U)[0], &block_size, &(*P)[0], &TMP[0],
	        &block_size, &info);

	for (int j = 0; j < block_size; j++) {
		for (int i = 0; i < block_size; i++) {
			(*L)[j + block_size * i] = TMP[i + block_size * j];
		}
	}
}
void RBMatrix::DTSTRF(blk_ptr A, blk_ptr L, blk_ptr U, piv_ptr P)
{
	P = rcp(new valarray<int>(2 * block_size + 1));
	valarray<double> tall(2 * block_size * block_size);
	tall[gslice(0, {(size_t) block_size, (size_t) block_size}, {1, (size_t) block_size})] = *U;
	tall[gslice(block_size * block_size, {(size_t) block_size, (size_t) block_size},
	            {1, (size_t) block_size})]
	= *A;
	int info;
	dgetrf_(&block_size, &block_size, &tall[0], &block_size, &(*P)[0], &info);
	*U = tall[gslice(0, {(size_t) block_size, (size_t) block_size}, {1, (size_t) block_size})];
	*L = tall[gslice(block_size * block_size, {(size_t) block_size, (size_t) block_size},
	                 {1, (size_t) block_size})];
}
void RBMatrix::DSSSSM(blk_ptr U, blk_ptr A, blk_ptr L, piv_ptr P){}
blk_ptr RBMatrix::blkCopy(Block in)
{
	blk_ptr C = rcp(new valarray<double>(in.block->size()));
	for (int i = 0; i < block_size; i++) {
        int in_i = i;
        if (in.flip_i){
			in_i = block_size - 1 - i;
		}
		for (int j = 0; j < block_size; j++) {
			int in_j = j;
			if (in.flip_j) {
				in_j = block_size - 1 - j;
			}
			(*C)[i + block_size * j] = (*in.block)[in_i + block_size * in_j];
		}
	}
	return C;
}
blk_ptr RBMatrix::getBlock(int i, int j){
	blk_ptr C;
	try{
		int local_i = range_map.at(i);
        C = block_cols.at(j/block_size).at(local_i).block;
	} catch (const out_of_range &oor) {
		// do nothing
    }
    return C;
}
void RBMatrix::lu(RCP<RBMatrix> &L, RCP<RBMatrix> &U)
{
	L = rcp(new RBMatrix(domain, block_size, num_blocks));
	U = rcp(new RBMatrix(domain, block_size, num_blocks));
    map<tuple<Block,Block>,blk_ptr> mm_comp;
	for (int k = 0; k < num_blocks; k++) {
		int range_k = range->getLocalElement(k * block_size);

		blk_ptr Akk = blkCopy(block_cols[k][range_k]);
		blk_ptr Lkk, Ukk;
		piv_ptr Pkk;
		DGETRF(Akk, Lkk, Ukk, Pkk);
		L->insertBlock(k * block_size, k * block_size, Lkk, false, false);
		U->insertBlock(k * block_size, k * block_size, Ukk, false, false);
		// get row of U's
		for (int j = k + 1; j < num_blocks; j++) {
			auto Akj_block = block_cols[j].find(range_k);
			if (Akj_block != block_cols[j].end()) {
				blk_ptr Ukj;
				blk_ptr Akj = blkCopy(Akj_block->second);
				DGESSM(Akj, Lkk, Ukj, Pkk);
				U->insertBlock(k * block_size, j * block_size, Ukj, false, false);
			}
		}
		// get Column of L's
		for (int i = k + 1; i < num_blocks; i++) {
			int  range_i   = range->getLocalElement(i * block_size);
			auto Aik_block = block_cols[k].find(range_i);
			if (Aik_block != block_cols[k].end()) {
				blk_ptr Lik;
				blk_ptr Aik = blkCopy(Aik_block->second);
				LSolve(Aik, Lik, Ukk, Pkk);
				L->insertBlock(i * block_size, k * block_size, Lik, false, false);
			}
		}
        //update rest of matrix
		for (int j = k + 1; j < num_blocks; j++) {
			int  range_k   = U->range_map[k * block_size];
			auto Ukj_block = U->block_cols[j].find(range_k);
			if (Ukj_block != U->block_cols[j].end()) {
				blk_ptr Ukj = Ukj_block->second.block;
				for (int i = k + 1; i < num_blocks; i++) {
					int  range_i   = L->range_map[i * block_size];
					auto Lik_block = L->block_cols[k].find(range_i);
					if (Lik_block != L->block_cols[k].end()) {
						blk_ptr Lik   = Lik_block->second.block;
						blk_ptr C;
						C             = rcp(new valarray<double>(block_size * block_size));
						char    trans = 'N';
						double  one   = -1;
						double  zero  = 0;
						dgemm_(&trans, &trans, &block_size, &block_size, &block_size, &one,
						       &(*Lik)[0], &block_size, &(*Ukj)[0], &block_size, &zero, &(*C)[0],
						       &block_size);
						insertBlock(i * block_size, j * block_size, C, false, false);
					}
				}
			}
		}
	}
}
void RBMatrix::ilu(RCP<RBMatrix> &L, RCP<RBMatrix> &U)
{
	L = rcp(new RBMatrix(domain, block_size, num_blocks));
	U = rcp(new RBMatrix(domain, block_size, num_blocks));
    map<tuple<Block,Block>,blk_ptr> mm_comp;
	for (int k = 0; k < num_blocks; k++) {
		int range_k = range->getLocalElement(k * block_size);

		blk_ptr Akk = blkCopy(block_cols[k][range_k]);
		blk_ptr Lkk, Ukk;
		piv_ptr Pkk;
		DGETRF(Akk, Lkk, Ukk, Pkk);
		L->insertBlock(k * block_size, k * block_size, Lkk, false, false);
		U->insertBlock(k * block_size, k * block_size, Ukk, false, false);
		// get row of U's
		for (int j = k + 1; j < num_blocks; j++) {
			auto Akj_block = block_cols[j].find(range_k);
			if (Akj_block != block_cols[j].end()) {
				blk_ptr Ukj;
				blk_ptr Akj = blkCopy(Akj_block->second);
				DGESSM(Akj, Lkk, Ukj, Pkk);
				U->insertBlock(k * block_size, j * block_size, Ukj, false, false);
			}
		}
		// get Column of L's
		for (int i = k + 1; i < num_blocks; i++) {
			int  range_i   = range->getLocalElement(i * block_size);
			auto Aik_block = block_cols[k].find(range_i);
			if (Aik_block != block_cols[k].end()) {
				blk_ptr Lik;
				blk_ptr Aik = blkCopy(Aik_block->second);
				LSolve(Aik, Lik, Ukk, Pkk);
				L->insertBlock(i * block_size, k * block_size, Lik, false, false);
			}
		}
        //update rest of matrix
		for (int j = k + 1; j < num_blocks; j++) {
			int  range_k   = U->range_map[k * block_size];
			blk_ptr Ukj       = U->getBlock(range_k, j * block_size);
			if (!Ukj.is_null()){
				for (int i = k + 1; i < num_blocks; i++) {
					int     range_i = L->range_map[i * block_size];
					blk_ptr Lik     = L->getBlock(i * block_size, range_i);
					if (!Lik.is_null() && getBlock(i * block_size, j * block_size).is_null()) {
						blk_ptr C;
						C             = rcp(new valarray<double>(block_size * block_size));
						char    trans = 'N';
						double  one   = -1;
						double  zero  = 0;
						dgemm_(&trans, &trans, &block_size, &block_size, &block_size, &one,
						       &(*Lik)[0], &block_size, &(*Ukj)[0], &block_size, &zero, &(*C)[0],
						       &block_size);
						insertBlock(i * block_size, j * block_size, C, false, false);
					}
				}
			}
		}
	}
}
ostream& operator<<(ostream& os, const RBMatrix& A) {
    os << scientific;
    int size= A.num_blocks*A.block_size;
    os.precision(15);
	os << "%%MatrixMarket matrix coordinate real general\n";
	os << size << ' ' << size << ' ' << A.nz << '\n';
	for (size_t index = 0; index < A.block_cols.size(); index++) {
		int start_j = A.domain->getGlobalElement(index * A.block_size);
		for (auto &p : A.block_cols[index]) {
			Block                 curr_block = p.second;
			std::valarray<double> curr_blk   = *curr_block.block;
			int                   start_i    = A.range->getGlobalElement(p.first);
			for (int iii = 0; iii < A.block_size; iii++) {
				int i = start_i + iii;
                int block_i = iii;
				if (curr_block.flip_i) {
					block_i = A.block_size - 1 - iii;
				}
				for (int jjj = 0; jjj < A.block_size; jjj++) {
					int j       = start_j + jjj;
					int block_j = jjj;
					if (curr_block.flip_j) {
						block_j = A.block_size - 1 - jjj;
					}
					os << i + 1 << ' ' << j + 1 << ' ' << curr_blk[block_i + A.block_size * block_j]
					   << "\n";
				}
			}
		}
	}
    return os;
}
