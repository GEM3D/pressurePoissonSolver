#ifndef FUNCTIONWRAPPER_H
#define FUNCTIONWRAPPER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
class FuncWrap // : virtual Belos::Operator
{
	public:
	Teuchos::RCP<vector_type> b;
	DomainCollection *        dc;
	FuncWrap(Teuchos::RCP<vector_type> b, DomainCollection *dc)
	{
		this->b  = b;
		this->dc = dc;
	}
	void Apply(const vector_type &x, vector_type &y) const
	{
		dc->solveWithInterface(x, y);
		y.update(1, *b, -1);
	}
};
namespace Belos
{
template <>
void OperatorTraits<double, vector_type, FuncWrap>::Apply(const FuncWrap &wrapper,
                                                          const vector_type &x, vector_type &y,
                                                          Belos::ETrans trans)
{
	wrapper.Apply(x, y);
};
}
#endif
