#include "DomainCollection.h"
#include <functional>
class Init
{
	public:
	static void initNeumann(DomainCollection<3> &dc, int n, Vec f, Vec exact,
	                        std::function<double(double, double, double)> ffun,
	                        std::function<double(double, double, double)> efun,
	                        std::function<double(double, double, double)> nfunx,
	                        std::function<double(double, double, double)> nfuny,
	                        std::function<double(double, double, double)> nfunz);
	static void initDirichlet(DomainCollection<3> &dc, int n, Vec f, Vec exact,
	                          std::function<double(double, double, double)> ffun,
	                          std::function<double(double, double, double)> efun);
};
