#include "DomainCollection.h"
#include <functional>
class Init
{
	public:
	static void initNeumann(DomainCollection &dc, int n, Vec f, Vec exact,
	                        std::function<double(double, double)> ffun,
	                        std::function<double(double, double)> efun,
	                        std::function<double(double, double)> nfunx,
	                        std::function<double(double, double)> nfuny);
	static void initDirichlet(DomainCollection &dc, int n, Vec f, Vec exact,
	                          std::function<double(double, double)> ffun,
	                          std::function<double(double, double)> efun);
	static void fillSolution(DomainCollection &dc, vector_type &u,
	                         std::function<double(double, double, double)> fun, double time);
};
