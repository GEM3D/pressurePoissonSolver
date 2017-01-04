#ifndef FUNCTIONWRAPPER_H
#define FUNCTIONWRAPPER_H
#include "Domain.h"
#include "MyTypeDefs.h"
#include <BelosEpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>
#include <Teuchos_RCP.hpp>
class FuncWrap // : virtual Belos::Operator
{
	public:
	Teuchos::RCP<vector_type> b;
	Domain *                  d;
	FuncWrap(Teuchos::RCP<vector_type> b, Domain *d)
	{
		this->b = b;
		this->d = d;
	}
	void Apply(const vector_type &x, vector_type &y) const
	{
		d->solveWithInterface(x, y);
		y.Update(1, *b, -1);
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
