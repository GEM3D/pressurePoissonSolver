#ifndef GMGHELPER_H
#define GMGHELPER_H
#include "DomainCollection.h"
#include "MatrixHelper.h"
#include "SchurHelper.h"
#include "GMGRestrictor.h"
#include "GMGInterpolator.h"
#include <petscpc.h>
class GMGHelper
{
	private:
	int                           num_levels;
	int                           top_level;
	std::vector<DomainCollection> levels;
	std::vector<SchurHelper>      shs;
	std::vector<PW<Vec>>          u_vectors;
	std::vector<PW<Vec>>          f_vectors;
	std::vector<PW<Vec>>          r_vectors;
    std::vector<GMGRestrictor*>  restrictors;
    std::vector<GMGInterpolator*>  interpolators;

	void apply(Vec f, Vec u);
	public:
	static int multiply(PC A, Vec f, Vec u)
	{
		GMGHelper *gh = nullptr;
		PCShellGetContext(A, (void **) &gh);
		VecScale(u, 0);
		gh->apply(f, u);
		return 0;
	}

	GMGHelper(int n, OctTree t, DomainCollection &dc,SchurHelper &sh);

	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};
#endif
