#include "Factory.h"
#include <MueLu.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_TpetraOperator_def.hpp>
Teuchos::RCP<op_type> Factory::getAmgPreconditioner(Teuchos::RCP<op_type>     A,
                                                    Teuchos::RCP<vector_type> xy)
{
	return MueLu::CreateTpetraPreconditioner(A, "mueluOptions.xml", xy);
}
Teuchos::RCP<matrix_type> Factory::blockInverse(Teuchos::RCP<matrix_type> A, int n) {}
