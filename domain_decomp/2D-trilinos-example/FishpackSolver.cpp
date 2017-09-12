#include "FishpackSolver.h"
extern "C" {
void hstcrt_(double *a, double *b, int *m, int *mbdcnd, double *bda, double *bdb, double *c,
             double *d, int *n, int *nbdcnd, double *bdc, double *bdd, double *elmbda, double *f,
             int *idimf, double *pertrb, int *ierror, double *w);
}
using namespace std;
FishpackSolver::FishpackSolver(Domain *dom)
{
	this->d = dom;
}
void FishpackSolver::solve()
{
	double            a = d->x_start;
	double            b = d->x_start + d->x_length;
	int               m = d->n;
	int               mbcdnd = -1;
	if (d->neumann && !(d->hasNbr(Side::east) || d->hasNbr(Side::west))) {
        mbcdnd = 3;
	} else if (d->neumann && !d->hasNbr(Side::west)){
        mbcdnd = 4;
	} else if (d->neumann && !d->hasNbr(Side::east)) {
        mbcdnd = 2;
	} else {
		mbcdnd = 1;
	}
    valarray<double> bda = d->boundary_west;
    valarray<double> bdb = d->boundary_east;

	double            c = d->y_start;
	double            d2 = d->y_start + d->y_length;
	int               n = d->n;
	int               nbcdnd = -1;
	if (d->neumann && !(d->hasNbr(Side::south) || d->hasNbr(Side::north))) {
        nbcdnd = 3;
	} else if (d->neumann && !d->hasNbr(Side::south)) {
        nbcdnd = 4;
	} else if (d->neumann && !d->hasNbr(Side::north)){
        nbcdnd = 2;
	} else {
		nbcdnd = 1;
	}
    valarray<double> bdc = d->boundary_south;
    valarray<double> bdd = d->boundary_north;

	double elmbda = 0;
	valarray<double> f      = d->f;
	int               idimf  = n;
    double pertrb=0;
    int ierror=0;
	valarray<double> w(13 * m + 4 * n + m * log2(n));
    w[0] = 30.5;
	hstcrt_(&a, &b, &m, &mbcdnd, &bda[0], &bdb[0], &c, &d2, &n, &nbcdnd, &bdc[0], &bdd[0], &elmbda,
	        &f[0], &idimf, &pertrb, &ierror, &w[0]);
	//cerr << ierror << ", " << pertrb << ", " << elmbda << endl;
	d->u=f;
}
