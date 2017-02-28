#ifndef REPEATBLOCKMATRIX_H
#define REPEATBLOCKMATRIX_H
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
#include <valarray>
#include <map>
#include <iostream>
class Block{
	public:
	bool flip_i;
	bool flip_j;

	Teuchos::RCP<std::valarray<double>> block;

    Block(){}
	Block(Teuchos::RCP<std::valarray<double>> block, bool flip_i, bool flip_j)
	{
		this->block  = block;
		this->flip_i = flip_i;
		this->flip_j = flip_j;
	}
	friend bool operator==(const Block &l, const Block &r)
	{
		return std::tie(l.flip_i, l.flip_j, l.block) == std::tie(r.flip_i, r.flip_j, r.block);
	}
	friend bool operator<(const Block &l, const Block &r)
	{
        double* left_ptr = &(*l.block)[0];
        double* right_ptr = &(*r.block)[0];
		return std::tie(l.flip_i, l.flip_j, left_ptr)
		       < std::tie(r.flip_i, r.flip_j, right_ptr);
	}
};
class RBMatrix : public Tpetra::Operator<>
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

	int block_size;
	int num_blocks;
    int nz = 0;

	public:
	int              skip_index   = -1;
	RBMatrix(Teuchos::RCP<map_type> map, int block_size, int num_blocks);

	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           double alpha = Teuchos::ScalarTraits<double>::one(),
	           double beta  = Teuchos::ScalarTraits<double>::zero()) const;
	void insertBlock(int i, int j, Teuchos::RCP<std::valarray<double>> block, bool flip_i,
	                 bool flip_j);
	void                         createRangeMap();
	Teuchos::RCP<const map_type> getDomainMap() const { return domain; }
	Teuchos::RCP<const map_type> getRangeMap() const { return range; }
	Teuchos::RCP<RBMatrix>       invBlockDiag();
	Teuchos::RCP<RBMatrix>       lu();
	void DGETRF(Teuchos::RCP<std::valarray<double>> A, Teuchos::RCP<std::valarray<double>> L,
	            Teuchos::RCP<std::valarray<double>> U, Teuchos::RCP<std::valarray<int>> P);
	void DGESSM(Teuchos::RCP<std::valarray<double>> A, Teuchos::RCP<std::valarray<double>> L,
	            Teuchos::RCP<std::valarray<int>> P, Teuchos::RCP<std::valarray<double>> U);
	void DTSTRF(Teuchos::RCP<std::valarray<double>> U, Teuchos::RCP<std::valarray<double>> A,
	            Teuchos::RCP<std::valarray<int>> P);
	void DSSSSM(Teuchos::RCP<std::valarray<double>> U, Teuchos::RCP<std::valarray<double>> A,
	            Teuchos::RCP<std::valarray<double>> L, Teuchos::RCP<std::valarray<int>> P);
	friend std::ostream &operator<<(std::ostream &os, const RBMatrix &A);
};
#endif
