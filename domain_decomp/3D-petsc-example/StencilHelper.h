#ifndef STENCILHELPER_H
#define STENCILHELPER_H
#include "Domain.h"
#include <valarray>
class StencilHelper
{
	public:
	virtual ~StencilHelper() {}
	virtual int row(int xi, int yi)        = 0;
	virtual int size(int xi, int yi)       = 0;
	virtual double *coeffs(int xi, int yi) = 0;
	virtual int *cols(int xi, int yi)      = 0;
};
class DirichletSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    stridex;
	int    stridey;

	public:
	DirichletSH(Domain &d, Side s)
	{
		double h   = 0;
		int    idx = d.id_global * d.n * d.n * d.n;
		switch (s) {
			case Side::west:
				h       = d.x_length / d.n;
				start   = idx;
				stridex = d.n;
				stridey = d.n * d.n;
				break;
			case Side::east:
				h       = d.x_length / d.n;
				start   = idx + (d.n - 1);
				stridex = d.n;
				stridey = d.n * d.n;
				break;
			case Side::south:
				h       = d.y_length / d.n;
				start   = idx;
				stridex = 1;
				stridey = d.n * d.n;
				break;
			case Side::north:
				h       = d.y_length / d.n;
				start   = idx + (d.n - 1) * d.n;
				stridex = 1;
				stridey = d.n * d.n;
				break;
			case Side::bottom:
				h       = d.z_length / d.n;
				start   = idx;
				stridex = 1;
				stridey = d.n;
				break;
			case Side::top:
				h       = d.z_length / d.n;
				start   = idx + (d.n - 1) * d.n * d.n;
				stridex = 1;
				stridey = d.n;
				break;
		}
		coeff = -1.0 / (h * h);
	}
	int row(int xi, int yi) { return start + stridex * xi + stridey * yi; }
	int size(int xi, int yi) { return 1; }
	double *coeffs(int xi, int yi) { return &coeff; }
	int *cols(int xi, int yi)
	{
		col = start + stridex * xi + stridey * yi;
		return &col;
	}
};
class NeumannSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    stridex;
	int    stridey;

	public:
	NeumannSH(Domain &d, Side s)
	{
		double h   = 0;
		int    idx = d.id_global * d.n * d.n * d.n;
		switch (s) {
			case Side::west:
				h       = d.x_length / d.n;
				start   = idx;
				stridex = d.n;
				stridey = d.n * d.n;
				break;
			case Side::east:
				h       = d.x_length / d.n;
				start   = idx + (d.n - 1);
				stridex = d.n;
				stridey = d.n * d.n;
				break;
			case Side::south:
				h       = d.y_length / d.n;
				start   = idx;
				stridex = 1;
				stridey = d.n * d.n;
				break;
			case Side::north:
				h       = d.y_length / d.n;
				start   = idx + (d.n - 1) * d.n;
				stridex = 1;
				stridey = d.n * d.n;
				break;
			case Side::bottom:
				h       = d.z_length / d.n;
				start   = idx;
				stridex = 1;
				stridey = d.n;
				break;
			case Side::top:
				h       = d.z_length / d.n;
				start   = idx + (d.n - 1) * d.n * d.n;
				stridex = 1;
				stridey = d.n;
				break;
		}
		coeff = 1.0 / (h * h);
	}
	int row(int xi, int yi) { return start + stridex * xi + stridey * yi; }
	int size(int xi, int yi) { return 1; }
	double *coeffs(int xi, int yi) { return &coeff; }
	int *cols(int xi, int yi)
	{
		col = start + stridex * xi + stridey * yi;
		return &col;
	}
};
class NormalSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    nbr_start;
	int    stridex;
	int    stridey;

	public:
	NormalSH(Domain &d, Side s)
	{
		double h       = 0;
		int    idx     = d.id_global * d.n * d.n * d.n;
		int    nbr_idx = d.globalNbr(s) * d.n * d.n * d.n;
		switch (s) {
			case Side::west:
				h         = d.x_length / d.n;
				start     = idx;
				nbr_start = nbr_idx + (d.n - 1);
				stridex   = d.n;
				stridey   = d.n * d.n;
				break;
			case Side::east:
				h         = d.x_length / d.n;
				start     = idx + (d.n - 1);
				nbr_start = nbr_idx;
				stridex   = d.n;
				stridey   = d.n * d.n;
				break;
			case Side::south:
				h         = d.y_length / d.n;
				start     = idx;
				nbr_start = nbr_idx + (d.n - 1) * d.n;
				stridex   = 1;
				stridey   = d.n * d.n;
				break;
			case Side::north:
				h         = d.y_length / d.n;
				start     = idx + (d.n - 1) * d.n;
				nbr_start = nbr_idx;
				stridex   = 1;
				stridey   = d.n * d.n;
				break;
			case Side::bottom:
				h         = d.z_length / d.n;
				start     = idx;
				nbr_start = nbr_idx + (d.n - 1) * d.n * d.n;
				stridex   = 1;
				stridey   = d.n;
				break;
			case Side::top:
				h         = d.z_length / d.n;
				start     = idx + (d.n - 1) * d.n * d.n;
				nbr_start = nbr_idx;
				stridex   = 1;
				stridey   = d.n;
				break;
		}
		coeff = 1.0 / (h * h);
	}
	int row(int xi, int yi) { return start + stridex * xi + stridey * yi; }
	int size(int xi, int yi) { return 1; }
	double *coeffs(int xi, int yi) { return &coeff; }
	int *cols(int xi, int yi)
	{
		col = nbr_start + stridex * xi + stridey * yi;
		return &col;
	}
};
StencilHelper *getStencilHelper(Domain &d, Side s)
{
	StencilHelper *retval = nullptr;
	if (d.hasNbr(s)) {
		retval = new NormalSH(d, s);
	} else {
		if (d.isNeumann(s)) {
						retval = new NeumannSH(d, s);
		} else {
			retval = new DirichletSH(d, s);
		}
	}
	return retval;
}
#endif
