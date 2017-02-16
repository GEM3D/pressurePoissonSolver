#include "RBMatrix.h"

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
RBMatrix::RBMatrix(Teuchos::RCP<map_type> map, int block_size, int num_blocks)
{
	domain = map;
	range  = map;

	this->num_blocks = num_blocks;
	this->block_size = block_size;

	block_cols = vector<std::map<int, Block>>(num_blocks);
}
void RBMatrix::apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode, double alpha,
                     double beta) const
{
	auto x_view = x.getLocalView<Kokkos::HostSpace>();
    y.putScalar(beta);
	auto y_view = y.getLocalView<Kokkos::HostSpace>();
	// column loop
	for (int index = 0; index < num_blocks; index++) {
		int              start_j = index * block_size;
		valarray<double> slot(block_size);
		valarray<double> slot_rev(block_size);
		for (int i = 0; i < block_size; i++) {
			slot[i] = x_view(start_j + i, 0);
		}
		for (int i = 0; i < block_size; i++) {
			slot_rev[i] = slot[block_size - 1 - i];
		}
		// go down the column
		for (auto &p : block_cols[index]) {
			int   start_i    = p.first;
			Block curr_block = p.second;
			for (int iii = 0; iii < block_size; iii++) {
				int block_i = iii;
				if (curr_block.flip_i) {
					block_i = block_size - 1 - iii;
				}
				valarray<double> &curr_blk  = *curr_block.block;
				valarray<double>  blk_slice = curr_blk[slice(block_i * block_size, block_size, 1)];

				int matrix_i = start_i + iii;
				if (curr_block.flip_j) {
					y_view(matrix_i, 0) += (slot_rev * blk_slice).sum();
				} else {
					y_view(matrix_i, 0) += (slot * blk_slice).sum();
				}
			}

		}
	}
}
void RBMatrix::insertBlock(int i, int j, RCP<valarray<double>> block, bool flip_i, bool flip_j)
{
	Block b(block, flip_i, flip_j);

    int index = j/block_size;
	if (block_cols[index].count(i) > 0) {
		Block first_block = block_cols[index][i];

		pair<Block, Block> bpair(first_block, b);
		pair<Block, Block> bpair_rev(b, first_block);

		if (combos.count(bpair) > 0) {
			block_cols[index][i] = combos[bpair];
		} else if (combos.count(bpair_rev) > 0) {
			block_cols[index][i] = combos[bpair_rev];
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
					new_blk[i * block_size + j] = first_blk[block_i * block_size + block_j];
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
					new_blk[i * block_size + j] += second_blk[block_i * block_size + block_j];
				}
			}
			block_cols[index][i] = new_block;
			combos[bpair]        = new_block;
		}
	} else {
        nz+=block->size();
		block_cols[index][i] = b;
	}
}
ostream& operator<<(ostream& os, const RBMatrix& A) {
    os << scientific;
    int size= A.num_blocks*A.block_size;
    os.precision(15);
	os << "%%MatrixMarket matrix coordinate real general\n";
	os << size << ' ' << size << ' ' << A.nz << '\n';
	for (size_t index = 0; index < A.block_cols.size(); index++) {
		int start_j = index * A.block_size;
		for (auto &p : A.block_cols[index]) {
			Block                 curr_block = p.second;
			std::valarray<double> curr_blk   = *curr_block.block;
			int                   start_i    = p.first;
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
					os << i + 1 << ' ' << j + 1 << ' ' << curr_blk[block_i * A.block_size + block_j]
					   << "\n";
				}
			}
		}
	}
    return os;
}
