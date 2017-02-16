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
	return;
}
void RBMatrix::insertBlock(int i, int j, RCP<valarray<double>> block, bool flip_i, bool flip_j)
{
	cerr << "\n inserting a block into " << i << "," << j << "\n";
	Block b(block, flip_i, flip_j);

    int index = j/block_size;
    cerr << "index is " << index << "\n";
	block_cols[index][i] = b;
}
