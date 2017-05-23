#ifndef FACTORY_H
#define FACTORY_H
#include "MyTypeDefs.h"
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosLSQRSolMgr.hpp>
#include <BelosLSQRSolMgr.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_RCP.hpp>
typedef Teuchos::RCP<Belos::LinearProblem<scalar_type, vector_type, Tpetra::Operator<scalar_type>>>
problem_type;
typedef Teuchos::RCP<Belos::SolverManager<scalar_type, vector_type, Tpetra::Operator<scalar_type>>>
solver_type;
typedef Tpetra::Operator<scalar_type> op_type;
class Factory
{
	public:
	static Teuchos::RCP<op_type> getAmgPreconditioner(Teuchos::RCP<op_type> A);
};
#endif
