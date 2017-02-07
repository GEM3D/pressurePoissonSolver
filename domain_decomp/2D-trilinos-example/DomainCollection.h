#include "Domain.h"
#include "MyTypeDefs.h"
#include "RBMatrix.h"
class DomainCollection
{
	private:
	std::map<int, Domain *> domains;
	Teuchos::RCP<const Teuchos::Comm<int>> comm;
	int                                    nx;
	int                                    ny;
	int                                    d_x;
	int                                    d_y;
	double                                 h_x;
	double                                 h_y;
	int                                    num_domains;
	int                                    num_global_domains;

	public:
	Teuchos::RCP<map_type>                 collection_map;
	Teuchos::RCP<map_type>                 collection_iface_map;
	Teuchos::RCP<map_type>                 matrix_map;
	Teuchos::RCP<map_type>                 iface_map;
	Teuchos::RCP<int_vector_type>          iface_info;

	DomainCollection(int low, int high, int nx, int ny, int d_x, int d_y, double h_x, double h_y,
	                 Teuchos::RCP<const Teuchos::Comm<int>> comm);

	void initDirichlet(std::function<double(double, double)> ffun,
	                   std::function<double(double, double)> gfun);

	void initNeumann(std::function<double(double, double)> ffun,
	                 std::function<double(double, double)> efun,
	                 std::function<double(double, double)> nfunx,
	                 std::function<double(double, double)> nfuny);

	void solveWithInterface(const vector_type &gamma, vector_type &diff);
	void                      generateMaps();
	void                      distributeIfaceInfo();
	double                    diffNorm();
	double diffNorm(double uavg, double eavg);
	double                    exactNorm();
	double exactNorm(double eavg);
	double                    uSum();
	double                    exactSum();
	Teuchos::RCP<matrix_type> formMatrix(Teuchos::RCP<map_type> map);
	Teuchos::RCP<matrix_type> formInvDiag(Teuchos::RCP<map_type> map);
	Teuchos::RCP<RBMatrix> formRBMatrix(Teuchos::RCP<map_type> map);
};
