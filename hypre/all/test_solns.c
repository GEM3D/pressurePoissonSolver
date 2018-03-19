#include "HYPRE_sstruct_ls.h"

#include "test_solns.h"

void qtrue_constant(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy)
{
    double c = 2;
    *q = c;
    *qx = 0;
    *qy = 0;
    *qxx = 0;
    *qyy = 0;
}

void qtrue_plane(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy)
{
    double a = 0;
    double b = 1;
    *q = a*x + b*y;
    *qx = a;
    *qy = b;
    *qxx = 0;
    *qyy = 0;
}

void qtrue_parabola(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy)
{
    double a = 1;
    double b = 0;
    *q = a*x*x + b*y*y;
    *qx = 2*a*x;
    *qy = 2*b*y;
    *qxx = 2*a;
    *qyy = 2*b;
}


void qtrue_trig(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy)
{
    double C = 2;
    double D = 2;
    double Cpi = C*M_PI;
    double Dpi = D*M_PI;
    *q = cos(Cpi*x)*cos(Dpi*y);
    *qx = -Cpi*sin(Cpi*x)*cos(Dpi*y);
    *qy = -Dpi*cos(Cpi*x)*sin(Dpi*y);
    *qxx = -Cpi*Cpi*(*q);
    *qyy = -Dpi*Dpi*(*q);
}
