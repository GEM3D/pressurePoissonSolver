#ifndef REPEATBLOCKMATRIX_H
#define REPEATBLOCKMATRIX_H
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
#include <iostream>
#include <map>
#include <valarray>
class Block
{
	public:
	bool   flip_i;
	bool   flip_j;
	double scale = 1;

	Teuchos::RCP<std::valarray<double>> block;

	Block() {}
	Block(Teuchos::RCP<std::valarray<double>> block, bool flip_i = false, bool flip_j = false,
	      double scale = 1)
	{
		this->block  = block;
		this->flip_i = flip_i;
		this->flip_j = flip_j;
		this->scale  = scale;
	}
	friend bool operator==(const Block &l, const Block &r)
	{
		return std::tie(l.flip_i, l.flip_j, l.block, l.scale)
		       == std::tie(r.flip_i, r.flip_j, r.block, r.scale);
	}
	friend bool operator<(const Block &l, const Block &r)
	{
		double *left_ptr  = &(*l.block)[0];
		double *right_ptr = &(*r.block)[0];
		return std::tie(l.flip_i, l.flip_j, left_ptr, l.scale)
		       < std::tie(r.flip_i, r.flip_j, right_ptr, r.scale);
	}
};
class RBMatrix : public Tpetra::Operator<scalar_type>
{
	private:
	Teuchos::RCP<map_type>         domain;
	Teuchos::RCP<map_type>         range;
	Teuchos::RCP<Tpetra::Export<>> exporter;
	std::map<int, int> range_map;
	int              curr_local_i = 0;
	int              local_skip_j = -1;
	int              local_skip_i = -1;
	std::vector<int> global_i;
	std::vector<std::map<int, Block>> block_cols;

	std::map<std::pair<Block, Block>, Block> combos;

	int  block_size;
	int  num_blocks;
	int  nz     = 0;
	bool zero_u = false;

	public:
	int skip_index = -1;

	RBMatrix(Teuchos::RCP<map_type> map, int block_size, int num_blocks);

	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
	           scalar_type beta  = Teuchos::ScalarTraits<scalar_type>::zero()) const;
	void insertBlock(int i, int j, Teuchos::RCP<std::valarray<double>> block, bool flip_i = false,
	                 bool flip_j = false, double scale = 1);
	void                                createRangeMap();
	Teuchos::RCP<const map_type>        getDomainMap() const { return domain; }
	Teuchos::RCP<const map_type>        getRangeMap() const { return range; }
	Teuchos::RCP<RBMatrix>              invBlockDiag();
	Teuchos::RCP<RBMatrix>              invBlockDiagR();
	Teuchos::RCP<std::valarray<double>> blkCopy(Block in);
	Teuchos::RCP<std::valarray<double>> getBlock(int i, int j);
	int                  getNumBlocks() { return num_blocks; }
	int                  getBlockSize() { return block_size; }
	friend std::ostream &operator<<(std::ostream &os, const RBMatrix &A);
};

#endif
