#ifndef GMGHELPER_H
#define GMGHELPER_H
#include "Cycle.h"
#include "DomainCollection.h"
#include "SchurHelper.h"
#include <petscpc.h>
namespace GMG
{
class Helper2d
{
	private:
	std::unique_ptr<Cycle> cycle;

	void apply(Vec f, Vec u);

	public:
	static int multiply(PC A, Vec f, Vec u)
	{
		Helper2d *gh = nullptr;
		PCShellGetContext(A, (void **) &gh);
		VecScale(u, 0);
		gh->apply(f, u);
		return 0;
	}

	Helper2d(int n, std::vector<std::shared_ptr<DomainCollection<2>>> domains,
	         std::shared_ptr<SchurHelper<2>> sh, std::string config_file);

	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};
} // namespace GMG
#endif
