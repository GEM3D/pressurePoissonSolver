#ifndef FUNCTIONWRAPPER_H
#define FUNCTIONWRAPPER_H
#include "SchurHelper.h"
#include <iostream>
class FuncWrap
{
	public:
	PW<Vec>           u;
	PW<Vec>           f;
	SchurHelper *     sh = nullptr;
	DomainCollection *dc = nullptr;
	FuncWrap()           = default;
	FuncWrap(SchurHelper *sh, DomainCollection *dc)
	{
		f        = dc->getNewDomainVec();
		u        = dc->getNewDomainVec();
		this->sh = sh;
		this->dc = dc;
	}
	static int multiply(Mat A, Vec x, Vec y)
	{
		FuncWrap *w = nullptr;
		MatShellGetContext(A, &w);
		w->sh->solveWithInterface(w->f, w->u, x, y);
		return 0;
	}
	PW_explicit<Mat> getMatrix()
	{
		PW<Mat> A;
		int     M = dc->num_global_interfaces * pow(dc->n, 2);
		int     m = dc->ifaces.size() * pow(dc->n, 2);
		MatCreateShell(MPI_COMM_WORLD, m, m, M, M, this, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) multiply);
		return A;
	}
};
#endif
