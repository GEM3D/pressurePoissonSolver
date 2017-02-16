#ifndef REPEATBLOCKMATRIX_H
#define REPEATBLOCKMATRIX_H
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
#include <valarray>
#include <map>
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
};
class RBMatrix : public Tpetra::Operator<>
{
	public:
	Teuchos::RCP<map_type> domain;
	Teuchos::RCP<map_type> range;
	std::vector<std::map<int, Block>> block_cols;

	int block_size;
	int num_blocks;

	RBMatrix(Teuchos::RCP<map_type> map,int block_size,int num_blocks);

	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           double alpha = Teuchos::ScalarTraits<double>::one(),
	           double beta  = Teuchos::ScalarTraits<double>::zero()) const;
	void insertBlock(int i, int j, Teuchos::RCP<std::valarray<double>> block, bool flip_i,
	                 bool flip_j);
	Teuchos::RCP<const map_type> getDomainMap() const { return domain; }
	Teuchos::RCP<const map_type> getRangeMap() const { return range; }
};
#endif
