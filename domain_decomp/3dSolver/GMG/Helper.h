#ifndef GMGHELPER_H
#define GMGHELPER_H
#include "Cycle.h"
#include "DomainCollection.h"
#include "SchurHelper.h"
#include "tpl/json.hpp"
#include <petscpc.h>
namespace GMG
{
class Helper
{
	private:
	std::unique_ptr<Cycle> cycle;

	void apply(Vec f, Vec u);

	public:
	static int multiply(PC A, Vec f, Vec u)
	{
		Helper *gh = nullptr;
		PCShellGetContext(A, (void **) &gh);
		VecScale(u, 0);
		gh->apply(f, u);
		return 0;
	}

	Helper(int n, std::vector<std::shared_ptr<DomainCollection<3>>> domains,
	       std::shared_ptr<SchurHelper<3>> sh, std::string config_file);

	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiply);
	}
};
} // namespace GMG
#endif
