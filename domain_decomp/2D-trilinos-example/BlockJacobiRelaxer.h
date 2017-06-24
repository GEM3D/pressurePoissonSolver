#ifndef BLOCKJACOBIRELAXER_H
#define BLOCKJACOBIRELAXER_H
#include "RBMatrix.h"
#include "OpShift.h"
#include "MyTypeDefs.h"
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
class BlockJacobiRelaxer : public Tpetra::Operator<scalar_type>
{
	public:
	Teuchos::RCP<RBMatrix>           D;
	Teuchos::RCP<Tpetra::Operator<>> A;
	int                              num_sweeps = 1;
	static double                    omega;
	BlockJacobiRelaxer(Teuchos::RCP<RBMatrix> A, Teuchos::RCP<single_vector_type> s = Teuchos::null)
	{
		D = A->invBlockDiag();
        if(s.is_null()){
			this->A = A;
		}else{
            Teuchos::RCP<OpShift> os = Teuchos::rcp(new OpShift(A, s));
			this->A         = os;
		}
	}
	void apply(const vector_type &x, vector_type &y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
	           scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
	           scalar_type beta  = Teuchos::ScalarTraits<scalar_type>::zero()) const
	{
		// D^-1(b-Rx);
		vector_type temp(x.getMap(),1);
		if (beta == 1.0) {
			for (int i = 0; i < num_sweeps; i++) {
				A->apply(y, temp);
				temp.update(1.0, x, -1.0);
				D->apply(temp, y, Teuchos::NO_TRANS, omega, 1.0+(1.0-omega));
			}
		} else {
			D->apply(x, y);
			for (int i = 0; i < num_sweeps - 1; i++) {
				A->apply(y, temp);
				temp.update(1.0, x, -1.0);
				D->apply(temp, y, Teuchos::NO_TRANS, omega, 1.0+(1.0-omega));
			}
		}
	}
	Teuchos::RCP<const map_type> getDomainMap() const { return D->getDomainMap(); }
	Teuchos::RCP<const map_type> getRangeMap() const { return D->getRangeMap(); }
};
double BlockJacobiRelaxer::omega = 1.0;
#endif

