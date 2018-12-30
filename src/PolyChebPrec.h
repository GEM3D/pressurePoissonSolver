#ifndef POLYCHEBPREC_H
#define POLYCHEBPREC_H
#include "DomainCollection.h"
#include "MatrixHelper.h"
#include "SchurHelper.h"
#include <petscpc.h>
class PolyChebPrec
{
	private:
	SchurHelper<3> *     sh;
	DomainCollection<3> *dc;
	double               interval = 0.95;
	std::vector<double>  coeffs
	= {4.472135954953655e+00, 5.675247900481234e+00, 3.601012922685066e+00, 2.284885928634731e+00,
	   1.449787551186771e+00, 9.199076055378766e-01, 5.836924189936992e-01, 3.703598469934007e-01,
	   2.349977690621489e-01, 1.491089055767314e-01, 9.461139059090561e-02, 6.003206306517687e-02,
	   3.809106471898141e-02, 2.416923786484517e-02, 1.533567161022980e-02, 1.628851184599676e-02};

	void apply(Vec f, Vec u);

	public:
	static int multiply(PC A, Vec b, Vec y)
	{
		PolyChebPrec *pc = nullptr;
		PCShellGetContext(A, (void **) &pc);
		pc->apply(b, y);
		return 0;
	}

	PolyChebPrec(SchurHelper<3> &sh, DomainCollection<3> &dc);

	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};
#endif
