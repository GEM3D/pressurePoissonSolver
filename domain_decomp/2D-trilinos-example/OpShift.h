#ifndef OPSHIFT_H
#define OPSHIFT_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
class OpShift : public Tpetra::Operator<scalar_type>
{
	public:
	Teuchos::RCP<matrix_type> A;
	Teuchos::RCP<single_vector_type> s;
    Teuchos::RCP<single_vector_type> ones;
	OpShift(Teuchos::RCP<matrix_type> A, Teuchos::RCP<single_vector_type> s)
	{
		this->A = A;
		this->s = s;
		ones    = Teuchos::rcp(new single_vector_type(s->getMap()));
		ones->putScalar(1);
        std::cerr << "Sum : "<<s->meanValue()*s->getGlobalLength()<<std::endl;
	}
	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
	           scalar_type beta  = Teuchos::ScalarTraits<scalar_type>::zero()) const
	{
        A->apply(x,y);
		Teuchos::ArrayView<double> dots(new double[s->getNumVectors()], s->getNumVectors());
		s->dot(y, dots);
		for (size_t i = 0; i < y.getNumVectors(); i++) {
            std::cerr << "Dot"<<i<<": " <<dots[i] << std::endl;
			y.getVectorNonConst(i)->update(dots[i], *ones, 1.0);
		}
	}
	Teuchos::RCP<const map_type> getDomainMap() const { return A->getDomainMap(); }
	Teuchos::RCP<const map_type> getRangeMap() const { return A->getRangeMap(); }
};
#endif
