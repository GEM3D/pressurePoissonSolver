#ifndef REPEATBLOCKMATRIX_H
#define REPEATBLOCKMATRIX_H
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
class RBMatrix : public Tpetra::Operator<>
{
	public:
	Teuchos::RCP<map_type> domain;
	Teuchos::RCP<map_type> range;
    std::vector<std::map<int,Teuchos::RCP<std::valarray<double>>>> block_cols;
    int block_size;
	RBMatrix(Teuchos::RCP<map_type> map,int block_size)
	{
        this->block_size = block_size;
	}
	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           double alpha = Teuchos::ScalarTraits<double>::one(),
	           double beta  = Teuchos::ScalarTraits<double>::zero()) const
	{

	}
    void insertBlock(int i,int j,Teuchos::RCP<std::valarray<double>> block){
    }
	Teuchos::RCP<const map_type> getDomainMap() const { return domain; }
	Teuchos::RCP<const map_type> getRangeMap() const { return range; }
};
#endif
