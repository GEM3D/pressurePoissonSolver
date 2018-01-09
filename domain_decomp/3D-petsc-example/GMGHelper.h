#ifndef GMGHELPER_H
#define GMGHELPER_H
#include "DomainCollection.h"
#include "FunctionWrapper.h"
#include "MatrixHelper.h"
class GMGHelper
{
	private:
	int                           num_levels;
	std::vector<DomainCollection> levels;
	std::vector<SchwarzPrec>      smoothers;
	std::vector<SchurHelper>      shs;
	std::vector<PW<Mat>>          matrices;
	std::vector<PW<Vec>>          u_vectors;
	std::vector<PW<Vec>>          f_vectors;
	std::vector<PW<Vec>>          r_vectors;

	void restrictForLevel(int level);
	void prolongateFromLevel(int level);
	void apply(Vec f, Vec u);
	static int multiply(PC A, Vec f, Vec u)
	{
		GMGHelper *gh = nullptr;
		PCShellGetContext(A, (void **) &gh);
		VecScale(u, 0);
		gh->apply(f, u);
		return 0;
	}

	public:
	GMGHelper(int n, OctTree t, std::shared_ptr<PatchSolver> solver,
	          std::shared_ptr<PatchOperator> op, std::shared_ptr<Interpolator> interpolator);

	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};
#endif
