#ifndef FUNCTIONWRAPPER_H
#define FUNCTIONWRAPPER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
class FuncWrap : public Tpetra::Operator<>
{
	public:
	Teuchos::RCP<vector_type> b;
	DomainCollection *        dc;
	FuncWrap(Teuchos::RCP<vector_type> b, DomainCollection *dc)
	{
		this->b  = b;
		this->dc = dc;
	}
	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           double alpha = Teuchos::ScalarTraits<double>::one(),
	           double beta  = Teuchos::ScalarTraits<double>::zero()) const
	{
		dc->solveWithInterface(x, y);
		y.update(1, *b, -1);

		auto y_view = y.getLocalView<Kokkos::HostSpace>();
		for (int i = 0; i < 400; i++) {
            std::cerr << y_view(i, 0) << " ";
		}
        std::cerr << std::endl;
	}
	Teuchos::RCP<const map_type> getDomainMap() const { return b->getMap(); }
	Teuchos::RCP<const map_type> getRangeMap() const { return b->getMap(); }
};
#endif
