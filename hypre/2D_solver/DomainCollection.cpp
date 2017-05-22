#include "DomainCollection.h"
#include <array>
#include <fstream>
#include <tuple>
using namespace std;

enum axis_enum { X_AXIS, Y_AXIS };
enum bc_enum { DIRICHLET, NEUMANN, REFINED };

DomainCollection::DomainCollection(DomainSignatureCollection dsc, int n)
{
	this->n            = n;
    this->dsc = dsc;
	num_global_domains = dsc.num_global_domains;
	for (auto p : dsc.domains) {
        //TODO change to shared ptr
		DomainSignature ds       = p.second;
		int             i        = ds.id;

		// create a domain
		domains[i]        =  Domain(ds, n);
	}
	// initialize hypre stuff
	HYPRE_SStructGridCreate(MPI_COMM_WORLD, 2, num_global_domains, &grid);
	for (int i = 0; i < num_global_domains; i++) {
		int low[2]  = {0, 0};
		int high[2] = {n-1, n-1};
		int ret = HYPRE_SStructGridSetExtents(grid, i, low, high);
		if (ret != 0) {
			throw std::runtime_error{"HYPRE_SStructGridSetExtents returned: " + ret};
        }
	}
    for(int i=0;i<num_global_domains;i++){
		int ret = HYPRE_SStructGridSetVariables(grid, i, 1, vartypes);
		if (ret != 0) {
			throw std::runtime_error{"HYPRE_SStructGridSetVariables returned: " + ret};
        }
	}
	for (auto &p : domains) {
		Domain &d        = p.second;
        d.setGridNbrs(grid);
    }
	HYPRE_SStructGridAssemble(grid);

	// define stencil
	int offsets[][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
	HYPRE_SStructStencilCreate(2, 5, &stencil_5pt);
	for (int entry = 0; entry < 5; ++entry) {
		HYPRE_SStructStencilSetEntry(stencil_5pt, entry, offsets[entry],0);
	}
}

void DomainCollection::initNeumann(function<double(double, double)> ffun,
                                   function<double(double, double)> efun,
                                   function<double(double, double)> nfunx,
                                   function<double(double, double)> nfuny, bool amr)
{
	neumann = true;
    this->amr = amr;
	for (auto &p : domains) {
		Domain &d        = p.second;

		// Generate RHS vector
		std::valarray<double> &f     = d.f;
		std::valarray<double> &exact = d.exact;

		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x           = d.x_start + d.h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + d.h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = efun(x, y);
				// west
				if (!d.hasNbr(Side::west)) {
					d.boundary_west[yi] = nfunx(d.x_start, y);
				}
				// east
				if (!d.hasNbr(Side::east)) {
					d.boundary_east[yi] = nfunx(d.x_start+d.x_length, y);
				}
				// south
				if (!d.hasNbr(Side::south)) {
					d.boundary_south[xi] = nfuny(x, d.y_start);
				}
				// north
				if (!d.hasNbr(Side::north)) {
					d.boundary_north[xi] = nfuny(x, d.y_start+d.y_length);
				}
			}
		}
		d.setNeumann();
	}
}

void DomainCollection::initDirichlet(function<double(double, double)> ffun,
                                     function<double(double, double)> gfun)
{
	for (auto &p : domains) {
		Domain &d        = p.second;

		// Generate RHS vector
		std::valarray<double> &f     = d.f;
		std::valarray<double> &exact = d.exact;
		// Use local indices to access the entries of f_data.
		for (int yi = 0; yi < n; yi++) {
			for (int xi = 0; xi < n; xi++) {
				double x           = d.x_start + d.h_x / 2.0 + d.x_length * xi / n;
				double y           = d.y_start + d.h_y / 2.0 + d.y_length * yi / n;
				f[yi * n + xi]     = ffun(x, y);
				exact[yi * n + xi] = gfun(x, y);
				if (!d.hasNbr(Side::north)) {
					d.boundary_north[xi] = gfun(x, d.y_start+d.y_length);
				}
				if (!d.hasNbr(Side::east)) {
					d.boundary_east[yi] = gfun(d.x_start+d.x_length, y);
				}
				if (!d.hasNbr(Side::south)) {
					d.boundary_south[xi] = gfun(x, d.y_start);
				}
				if (!d.hasNbr(Side::west)) {
					d.boundary_west[yi] = gfun(d.x_start, y);
				}
			}
		}
	}
}
void DomainCollection::formMatrix()
{
	// graph
	HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);
	if (use_parcsr) {
		HYPRE_SStructGraphSetObjectType(graph, HYPRE_PARCSR);
	}
	for (int i = 0; i < num_global_domains; i++) {
		HYPRE_SStructGraphSetStencil(graph, i, 0, stencil_5pt);
	}
    // add amr stencil entries
	for (auto &p : domains) {
		Domain &d        = p.second;
        d.setAmrStencil(graph);
    }
	HYPRE_SStructGraphAssemble(graph);


	HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
	if (use_parcsr) {
		HYPRE_SStructMatrixSetObjectType(A, HYPRE_PARCSR);
	}
	HYPRE_SStructMatrixInitialize(A);
	for (auto &p : domains) {
		Domain &d = p.second;
		d.setMatrixCoeffs(A);
	}
	HYPRE_SStructMatrixAssemble(A);
	if (use_parcsr) {
		HYPRE_SStructMatrixGetObject(A, (void **) &par_A);
	}
}
void DomainCollection::initVectors()
{
	HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &b);
	HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &x);
	if (use_parcsr) {
		HYPRE_SStructVectorSetObjectType(b, HYPRE_PARCSR);
		HYPRE_SStructVectorSetObjectType(x, HYPRE_PARCSR);
	}
	HYPRE_SStructVectorInitialize(b);
	HYPRE_SStructVectorInitialize(x);
	for (auto &p : domains) {
		Domain &d = p.second;
		d.fillRHS(b);
		//d.fillLHS(x);
	}
	HYPRE_SStructVectorAssemble(b);
	HYPRE_SStructVectorAssemble(x);
	if (use_parcsr) {
		HYPRE_SStructVectorGetObject(b, (void **) &par_b);
		HYPRE_SStructVectorGetObject(x, (void **) &par_x);
	}
}
void DomainCollection::saveResult()
{
	for (auto &p : domains) {
		Domain &d = p.second;
		//d.fillExact(x);
		// d.fillLHS(x);
	}
	HYPRE_SStructVectorAssemble(x);
	for (auto &p : domains) {
		Domain &d = p.second;
		d.saveLHS(x);
	}

	HYPRE_SStructMatrixMatvec(-1, A, x, 1, b);
	for (auto &p : domains) {
		Domain &d = p.second;
		d.saveResid(b);
	}
	HYPRE_SStructMatrixMatvec(1, A, x, 0, b);
	for (auto &p : domains) {
		Domain &d = p.second;
		d.saveAU(b);
	}
}
double DomainCollection::diffNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second.diffNorm(), 2);
	}
	return sqrt(result);
}
double DomainCollection::diffNorm(double uavg, double eavg)
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second.diffNorm(uavg, eavg), 2);
	}
	return sqrt(result);
}
double DomainCollection::exactNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second.exactNorm(), 2);
	}
	return sqrt(result);
}
double DomainCollection::fNorm()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second.fNorm(), 2);
	}
	return sqrt(result);
}
double DomainCollection::exactNorm(double eavg)
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second.exactNorm(eavg), 2);
	}
	return sqrt(result);
}
double DomainCollection::residual()
{
	double result = 0;
	for (auto &p : domains) {
		result += pow(p.second.residual(), 2);
	}
	return sqrt(result);
}
double DomainCollection::integrateF()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second.integrateF();
	}
	return sum;
}
double DomainCollection::integrateBoundaryFlux()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second.integrateBoundaryFlux();
	}
	return sum;
}
double DomainCollection::area()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second.area();
	}
	return sum;
}
double DomainCollection::integrateU()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second.integrateU();
	}
	return sum;
}
double DomainCollection::integrateExact()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second.integrateExact();
	}
	return sum;
}
double DomainCollection::integrateAU()
{
	double sum = 0;
	for (auto &p : domains) {
		sum += p.second.integrateAU();
	}
	return sum;
}

void DomainCollection::outputSolution(std::ostream &os)
{
	int num_i, num_j, d_x;
	if (amr) {
		num_i = n * sqrt(num_global_domains / 5);
		num_j = n * sqrt(num_global_domains / 5);
		d_x   = sqrt(num_global_domains / 5);
	} else {
		num_i = n * sqrt(num_global_domains);
		num_j = n * sqrt(num_global_domains);
		d_x   = sqrt(num_global_domains);
	}
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = domain_i * d_x + domain_j;
			os << domains[id].u[internal_i * n + internal_j] << '\n';
		}
	}
}
void DomainCollection::outputSolutionRefined(std::ostream &os)
{
	int num_i = 2 * n * sqrt(num_global_domains / 5);
	int num_j = 2 * n * sqrt(num_global_domains / 5);
	int d_x   = 2 * sqrt(num_global_domains / 5);
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = d_x * d_x / 4 + domain_i * d_x + domain_j;
			os << domains[id].u[internal_i * n + internal_j] << '\n';
		}
	}
}
void DomainCollection::outputResidual(std::ostream &os)
{
	int num_i, num_j, d_x;
	if (amr) {
		num_i = n * sqrt(num_global_domains / 5);
		num_j = n * sqrt(num_global_domains / 5);
		d_x   = sqrt(num_global_domains / 5);
	} else {
		num_i = n * sqrt(num_global_domains);
		num_j = n * sqrt(num_global_domains);
		d_x   = sqrt(num_global_domains);
	}
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = domain_i * d_x + domain_j;
            Domain &d = domains[id];
			os << d.resid[internal_i * n + internal_j]*d.h_x*d.h_y << '\n';
		}
	}
}
void DomainCollection::outputResidualRefined(std::ostream &os)
{
	int num_i = 2 * n * sqrt(num_global_domains / 5);
	int num_j = 2 * n * sqrt(num_global_domains / 5);
	int d_x   = 2 * sqrt(num_global_domains / 5);
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = d_x * d_x / 4 + domain_i * d_x + domain_j;
            Domain &d = domains[id];
			os << d.resid[internal_i * n + internal_j]*d.h_x*d.h_y << '\n';
		}
	}
}
void DomainCollection::outputError(std::ostream &os)
{
	int num_i, num_j, d_x;
	if (amr) {
		num_i = n * sqrt(num_global_domains / 5);
		num_j = n * sqrt(num_global_domains / 5);
		d_x   = sqrt(num_global_domains / 5);
	} else {
		num_i = n * sqrt(num_global_domains);
		num_j = n * sqrt(num_global_domains);
		d_x   = sqrt(num_global_domains);
	}
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = domain_i * d_x + domain_j;
			os << domains[id].error[internal_i * n + internal_j] << '\n';
		}
	}
}
void DomainCollection::outputErrorRefined(std::ostream &os)
{
	int num_i = 2 * n * sqrt(num_global_domains / 5);
	int num_j = 2 * n * sqrt(num_global_domains / 5);
	int d_x   = 2 * sqrt(num_global_domains / 5);
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = d_x * d_x / 4 + domain_i * d_x + domain_j;
			os << domains[id].error[internal_i * n + internal_j] << '\n';
		}
	}
}
void DomainCollection::outputClaw()
{
	ofstream     t_file("fort.t0000");
	const string tab = "\t";
	t_file << 0.0 << tab << "time" << endl;
	t_file << 3 << tab << "meqn" << endl;
	t_file << num_global_domains << tab << "ngrids" << endl;
	t_file << 2 << tab << "num_aux" << endl;
	t_file << 2 << tab << "num_dim" << endl;
	t_file.close();
	ofstream q_file("fort.q0000");

	q_file.precision(10);
	q_file << scientific;
	for (auto &p : domains) {
		Domain &d = p.second;
		d.outputClaw(q_file);
	}
	q_file.close();
}
