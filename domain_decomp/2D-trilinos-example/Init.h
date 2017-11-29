#include "DomainSignatureCollection.h"
#include <functional>
class Init
{
	public:
	static void initNeumann(DomainSignatureCollection &dsc, int n, double *f_vec, double *exact_vec,
	                        std::function<double(double, double)> ffun,
	                        std::function<double(double, double)> efun,
	                        std::function<double(double, double)> nfunx,
	                        std::function<double(double, double)> nfuny);
	static void initDirichlet(DomainSignatureCollection &dsc, int n, double *f_vec,
	                          double *exact_vec, std::function<double(double, double)> ffun,
	                          std::function<double(double, double)> efun);
};
