#include "RBMatrix.h"


extern "C" {

// LU decomoposition of a general matrix
void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *lda, int *IPIV, double *B, int *LDB,
             int *INFO);
void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA,
            double *B, int *LDB, double *BETA, double *C, int *LDC);
void dsyrk_(char*,char*,int*,int*,double*,double*,int*,double*,double*,int*);

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
void RBMatrix::apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode, scalar_type alpha,
                     scalar_type beta) const
{
    vector_type my_y(range, 1);
	y.scale(beta);
	my_y.doImport(y, *exporter, Tpetra::CombineMode::ADD);
	auto x_view = x.getLocalView<Kokkos::HostSpace>();
	auto y_view = my_y.getLocalView<Kokkos::HostSpace>();
	// column loop
	for (int index = 0; index < num_blocks; index++) {
		int              start_j = index * block_size;
		// go down the column
		for (auto &p : block_cols[index]) {
			int   start_i    = p.first;
			Block curr_block = p.second;
			for (int jjj = 0; jjj < block_size; jjj++) {
				int block_j = jjj;
				if (curr_block.flip_j) {
					block_j = block_size - 1 - jjj;
				}
				valarray<double> &curr_blk = *curr_block.block;

				if (curr_block.flip_i) {
					for (int i = 0; i < block_size; i++) {
						y_view(start_i + i, 0)
						+= x_view(start_j + jjj, 0)
						   * curr_blk[block_j * block_size + (block_size - 1 - i)];
					}
				} else {
					for (int i = 0; i < block_size; i++) {
						y_view(start_i + i, 0)
						+= x_view(start_j + jjj, 0) * curr_blk[block_j * block_size + i];
					}
				}
			}
		}
	}
    my_y.scale(alpha);
	y.doExport(my_y, *exporter, Tpetra::CombineMode::ADD);
}
void RBMatrix::insertBlock(int i, int j, RCP<valarray<double>> block, bool flip_i, bool flip_j,
                           double scale)
{
    i *=block_size;
    j *=block_size;
	int   local_j = domain->getLocalElement(j);
	int   local_i = -1;
	try{
		local_i = range_map.at(i);
	} catch (const out_of_range &oor) {
		local_i = curr_local_i;
		for (int x = 0; x < block_size; ++x) {
			global_i.push_back(i+x);
		}
		curr_local_i += block_size;
		range_map[i] = local_i;
	}

	Block b(block, flip_i, flip_j,scale);

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
					new_blk[i + block_size * j]
					= first_blk[block_i + block_size * block_j] * first_block.scale;
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
					new_blk[i + block_size * j]
					+= second_blk[block_i + block_size * block_j] * b.scale;
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
	exporter = rcp(new Tpetra::Export<>(range, domain));
}
RCP<RBMatrix> RBMatrix::invBlockDiag()
{
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
				try {
					blk_inv_ptr = computed.at(&(*curr_block.block)[0]);
				} catch (const out_of_range &oor) {
					// invert this block
					blk_inv_ptr               = blkCopy(curr_block);
					valarray<double> &blk_inv = *blk_inv_ptr;
                    for(int i=0;i<block_size;i++){
						blk_inv[i + i * block_size] += 2;
					}

					// compute inverse of block
					valarray<int>    ipiv(block_size + 1);
					int              lwork = block_size * block_size;
					valarray<double> work(lwork);
					int              info;
					dgetrf_(&block_size, &block_size, &blk_inv[0], &block_size, &ipiv[0], &info);
					dgetri_(&block_size, &blk_inv[0], &block_size, &ipiv[0], &work[0], &lwork,
					        &info);

					// insert the block
					computed[&(*curr_block.block)[0]] = blk_inv_ptr;
				}
				Inv->insertBlock(global_i / block_size, global_j / block_size, blk_inv_ptr);
			}
		}
	}
	Inv->createRangeMap();
	return Inv;
}
RCP<RBMatrix> RBMatrix::invBlockDiagR()
{
	RCP<RBMatrix>    R    = rcp(new RBMatrix(domain, block_size, num_blocks));
	for (int index = 0; index < num_blocks; index++) {
		int              start_j  = index * block_size;
		int              global_j = domain->getGlobalElement(start_j);
		// go down the column
		for (auto &p : block_cols[index]) {
			int   start_i    = p.first;
			int   global_i   = range->getGlobalElement(start_i);
			if (global_j != global_i) {
				Block                 curr = p.second;
				R->insertBlock(global_i / block_size, global_j / block_size, curr.block,
				               curr.flip_i, curr.flip_j);
			}
		}
	}
	R->createRangeMap();
	return R;
}
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
			(*C)[i + block_size * j] = (*in.block)[in_i + block_size * in_j] * in.scale;
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
ostream& operator<<(ostream& os, const RBMatrix& A) {
    os << scientific;
    int size= A.num_blocks*A.block_size;
    os.precision(15);
	os << "%%MatrixMarket matrix coordinate real general\n";
	os << size << ' ' << size << ' ' << A.nz << '\n';
	for (size_t index = 0; index < A.block_cols.size(); index++) {
		int start_j = A.domain->getGlobalElement(index * A.block_size);
        for(int local_i=0;local_i<A.block_size*A.num_blocks;local_i+=A.block_size){
			int start_i = A.range->getGlobalElement(local_i);
			try {
				Block                 curr_block = A.block_cols[index].at(local_i);
				std::valarray<double> curr_blk   = *curr_block.block;
				for (int iii = 0; iii < A.block_size; iii++) {
					int i       = start_i + iii;
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
						if (i == j) {
							os << i + 1 << ' ' << j + 1 << ' '
							   << (curr_blk[block_i + A.block_size * block_j])
							      * curr_block.scale
							   << "\n";
						} else {
							os << i + 1 << ' ' << j + 1 << ' '
							   << (curr_blk[block_i + A.block_size * block_j])
							      * curr_block.scale
							   << "\n";
						}
					}
				}
			} catch (out_of_range oor) {
			}
		}
	}
	return os;
}
