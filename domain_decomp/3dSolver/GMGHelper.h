#ifndef GMGHELPER_H
#define GMGHELPER_H
#include "DomainCollection.h"
#include "GMGInterpolator.h"
#include "GMGRestrictor.h"
#include "InterLevelComm.h"
#include "MatrixHelper.h"
#include "SchurHelper.h"
#include <petscpc.h>
class GMGHelper
{
	private:
	int                                            num_levels;
	int                                            top_level;
	std::vector<std::shared_ptr<DomainCollection>> levels;
	std::vector<SchurHelper>                       shs;
	std::vector<PW<Vec>>                           u_vectors;
	std::vector<PW<Vec>>                           f_vectors;
	std::vector<PW<Vec>>                           r_vectors;
	std::vector<PW<Mat>>                           mats;
	std::vector<std::shared_ptr<GMGRestrictor>>    restrictors;
	std::vector<std::shared_ptr<GMGInterpolator>>  interpolators;
	std::vector<std::shared_ptr<InterLevelComm>>   comms;

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

	GMGHelper(int n, OctTree t, std::shared_ptr<DomainCollection> dc, SchurHelper &sh);

	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};
#endif
