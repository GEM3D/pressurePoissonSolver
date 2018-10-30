#include "PolyChebPrec.h"
#include <iostream>
using namespace std;
PolyChebPrec::PolyChebPrec(SchurHelper &sh, DomainCollection<3> &dc)
{
	this->sh = &sh;
	this->dc = &dc;
}
void PolyChebPrec::apply(Vec b, Vec y)
{
	PW<Vec> f = dc->getNewDomainVec();
	PW<Vec> u = dc->getNewDomainVec();
	PW<Vec> bk;
	PW<Vec> bk1;
	PW<Vec> bk2;
	VecDuplicate(b, &bk);
	VecDuplicate(b, &bk1);
	VecDuplicate(b, &bk2);
	for (int i = coeffs.size() - 1; i > 0; i--) {
		sh->solveAndInterpolateWithInterface(f, u, bk1, bk);
		VecAXPBY(bk, -2, 4 / interval, bk1);
		VecAXPBYPCZ(bk, coeffs[i], -1.0, 1.0, b, bk2);
		PW<Vec> tmp = bk2;
		bk2         = bk1;
		bk1         = bk;
		bk          = tmp;
	}
	sh->solveAndInterpolateWithInterface(f, u, bk1, y);
	VecAXPBY(y, -1, 2 / interval, bk1);
	VecAXPBYPCZ(y, coeffs[0], -1.0, 1.0, b, bk2);
}
