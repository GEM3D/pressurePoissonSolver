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
				int block_i = iii;
				if (curr_block.flip_i) {
					block_i = block_size - 1 - iii;
				}
				valarray<double> &curr_blk  = *curr_block.block;

				int matrix_i = start_i + iii;
				if (curr_block.flip_j) {
					for (int i = 0; i < block_size; i++) {
						y_view(matrix_i, 0)
						+= x_view(start_j + i,0) * curr_blk[block_i * block_size+block_size-1 - i];
					}
				} else {
					for (int i = 0; i < block_size; i++) {
						y_view(matrix_i, 0)
						+= x_view(start_j + i,0) * curr_blk[block_i * block_size + i];
					}
				}
			}

		}
	}
	Tpetra::Export<> exporter(range, domain);
	y.doExport(my_y, exporter, Tpetra::CombineMode::ADD);
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
