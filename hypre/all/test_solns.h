#ifndef TEST_SOLNS_H
#define TEST_SOLNS_H


#include "HYPRE_sstruct_ls.h"   /* For M_PI */

#define  PI2 M_PI*M_PI


typedef void (*solution_t)(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy);


void qtrue_constant(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy);
void qtrue_plane(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy);
void qtrue_parabola(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy);
void qtrue_trig(double x, double y, double *q, double *qx, double *qy, double *qxx, double *qyy);

#endif
