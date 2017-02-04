#include "Domain.h"
#include "MyTypeDefs.h"
class DomainCollection
{
	private:
	std::map<int, Domain *> domains;
	Teuchos::RCP<map_type>                 collection_map;
	Teuchos::RCP<const Teuchos::Comm<int>> comm;
	int                                    nx;
	int                                    ny;
	double                                 h_x;
	double                                 h_y;
	int                                    num_domains;
	int                                    num_global_domains;

	public:
	DomainCollection(int low, int high, int nx, int ny, int d_x, int d_y, double h_x, double h_y,
	                 Teuchos::RCP<const Teuchos::Comm<int>> comm,
	                 std::function<double(double, double)> ffun,
	                 std::function<double(double, double)> gfun);
	void solveWithInterface(const vector_type &gamma, vector_type &diff);
    void generateMaps();
	double                    diffNorm();
	double                    exactNorm();
	Teuchos::RCP<matrix_type> formMatrix(Teuchos::RCP<map_type> map);
};
