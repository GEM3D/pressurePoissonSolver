#ifndef GMGSCHURHELPER_H
#define GMGSCHURHELPER_H
#include "DomainCollection.h"
#include "MatrixHelper.h"
#include "SchurHelper.h"
#include <petscpc.h>
class GMGSchurHelper
{
	private:
	int                           num_levels;
	int                           top_level;
	std::vector<DomainCollection> levels;
	std::vector<SchurHelper>      shs;
	std::vector<BlockJacobiSmoother>      jacobis;
	std::vector<PW<Vec>>          u_vectors;
	std::vector<PW<Vec>>          f_vectors;
	std::vector<PW<Vec>>          r_vectors;
	std::vector<PW<Vec>>          b_vectors;
	std::vector<PW<Vec>>          g_vectors;

	void restrictForLevel(int level);
	void prolongateFromLevel(int level);
	void apply(Vec f, Vec u);
	public:
	static int multiply(PC A, Vec f, Vec u)
	{
		GMGSchurHelper *gh = nullptr;
		PCShellGetContext(A, (void **) &gh);
		VecScale(u, 0);
		gh->apply(f, u);
		return 0;
	}

	GMGSchurHelper(int n, OctTree t, DomainCollection &dc,SchurHelper &sh);

	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};
#endif
