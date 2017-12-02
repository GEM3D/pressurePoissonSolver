#ifndef FUNCTIONWRAPPER_H
#define FUNCTIONWRAPPER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
class FuncWrap : public Tpetra::Operator<scalar_type>
{
	public:
	Teuchos::RCP<vector_type> b;
	DomainCollection *        dc;
	FuncWrap(Teuchos::RCP<vector_type> b, DomainCollection *dc)
	{
		this->b  = b;
		this->dc = dc;
	}
	void setB(Teuchos::RCP<vector_type> b) { this->b = b; }
	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
	           scalar_type beta  = Teuchos::ScalarTraits<scalar_type>::zero()) const
	{
		for (size_t i = 0; i < x.getNumVectors(); i++) {
			auto sx = *x.getVector(i);
			auto sy = *y.getVectorNonConst(i);
			dc->solveWithInterface(sx, sy);
			sy.update(1, *b, -1);
		}
	}
	Teuchos::RCP<const map_type> getDomainMap() const { return b->getMap(); }
	Teuchos::RCP<const map_type> getRangeMap() const { return b->getMap(); }
};
#endif
