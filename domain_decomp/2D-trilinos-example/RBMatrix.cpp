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
    cerr << "Block size is "<< block_size << "\n";
    cerr << "Num Blocks is "<< num_blocks << "\n";
}
void RBMatrix::apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode, double alpha,
                     double beta) const
{
    cerr << "apply has been called with alpha: " << alpha << " beta: " << beta << "\n";
	auto x_view = x.getLocalView<Kokkos::HostSpace>();
    y.putScalar(0);
	auto y_view = y.getLocalView<Kokkos::HostSpace>();
    //column loop
    for(size_t index = 0;index<block_cols.size();index++){
		int start_j = index * block_size;
		valarray<double> slot(block_size);
		valarray<double> slot_rev(block_size);
		for (int i = 0; i < block_size; i++) {
			slot[i]     = x_view(start_j + i, 0);
			slot_rev[i] = x_view(start_j + block_size - 1 - i, 0);
		}
		// go down the column
		for (auto &p : block_cols[index]) {
			Block             curr_block    = p.second;
			valarray<double> *curr_slot_ptr = &slot;
			if (curr_block.flip_j) {
				curr_slot_ptr = &slot_rev;
			}
			valarray<double> &curr_slot = *curr_slot_ptr;
			for (int index = 0; index < block_size; index++) {
				int block_i = index;
				if (curr_block.flip_i) {
					block_i = block_size - 1 - index;
				}
				valarray<double> &curr_blk = *curr_block.block;
				valarray<double> blk_slice = curr_blk[slice(block_i * block_size, block_size, 1)];

				int               matrix_i  = p.first + index;
				y_view(matrix_i, 0) += (curr_slot * blk_slice).sum();
			}
		}
	}
    cerr<<"norm of resulting vector " << y.getVector(0)->norm2()<<"\n";
}
void RBMatrix::insertBlock(int i, int j, RCP<valarray<double>> block, bool flip_i, bool flip_j)
{
	cerr << "inserting a block into " << i << "," << j << "\n";
	Block b(block, flip_i, flip_j);

    int index = j/block_size;
    cerr << "index is " << index << "\n\n";
	if (block_cols[index].count(i) > 0) {
		Block first_block = block_cols[index][i];

		pair<Block, Block> bpair(first_block, b);
		pair<Block, Block> bpair_rev(b, first_block);

		if (combos.count(bpair) > 0) {
			block_cols[index][i] = combos[bpair];
			cerr << "combo block reused \n\n";
		} else if (combos.count(bpair_rev) > 0) {
			block_cols[index][i] = combos[bpair_rev];
			cerr << "combo block reused \n\n";
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
			cerr << "new combo block created! \n\n";
		}
	} else {
		// do nothing
		block_cols[index][i] = b;
	}
}
