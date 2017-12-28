#ifndef STENCILHELPER_H
#define STENCILHELPER_H
#include "Domain.h"
#include <valarray>
class StencilHelper
{
	public:
	virtual ~StencilHelper() {}
	virtual int row(int i)        = 0;
	virtual int size(int i)       = 0;
	virtual double *coeffs(int i) = 0;
	virtual int *cols(int i)      = 0;
};
class DirichletSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    stride;

	public:
	DirichletSH(Domain &d, Side s)
	{
		double h   = 0;
		int    idx = d.id_global * d.n * d.n;
		switch (s) {
			case Side::north:
				h      = d.x_length / d.n;
				start  = idx + (d.n - 1) * d.n;
				stride = 1;
				break;
			case Side::east:
				h      = d.y_length / d.n;
				start  = idx + (d.n - 1);
				stride = d.n;
				break;
			case Side::south:
				h      = d.x_length / d.n;
				start  = idx;
				stride = 1;
				break;
			case Side::west:
				h      = d.y_length / d.n;
				start  = idx;
				stride = d.n;
				break;
            default:
                break;
		}
		coeff = -1.0 / (h * h);
	}
	int row(int i) { return start + stride * i; }
	int size(int i) { return 1; }
	double *coeffs(int i) { return &coeff; }
	int *cols(int i)
	{
		col = start + stride * i;
		return &col;
	}
};
class NeumannSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    stride;

	public:
	NeumannSH(Domain &d, Side s)
	{
		double h   = 0;
		int    idx = d.id_global * d.n * d.n;
		switch (s) {
			case Side::north:
				h      = d.x_length / d.n;
				start  = idx + (d.n - 1) * d.n;
				stride = 1;
				break;
			case Side::east:
				h      = d.y_length / d.n;
				start  = idx + (d.n - 1);
				stride = d.n;
				break;
			case Side::south:
				h      = d.x_length / d.n;
				start  = idx;
				stride = 1;
				break;
			case Side::west:
				h      = d.y_length / d.n;
				start  = idx;
				stride = d.n;
				break;
            default:
                break;
		}
		coeff = 1.0 / (h * h);
	}
	int row(int i) { return start + stride * i; }
	int size(int i) { return 1; }
	double *coeffs(int i) { return &coeff; }
	int *cols(int i)
	{
		col = start + stride * i;
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
	int    stride;

	public:
	NormalSH(Domain &d, Side s)
	{
		double h       = 0;
		int    idx     = d.id_global * d.n * d.n;
		int    nbr_idx = d.globalNbr(s) * d.n * d.n;
		switch (s) {
			case Side::north:
				h         = d.x_length / d.n;
				start     = idx + (d.n - 1) * d.n;
				nbr_start = nbr_idx;
				stride    = 1;
				break;
			case Side::east:
				h         = d.y_length / d.n;
				start     = idx + (d.n - 1);
				nbr_start = nbr_idx;
				stride    = d.n;
				break;
			case Side::south:
				h         = d.x_length / d.n;
				start     = idx;
				nbr_start = nbr_idx + (d.n - 1) * d.n;
				stride    = 1;
				break;
			case Side::west:
				h         = d.y_length / d.n;
				start     = idx;
				nbr_start = nbr_idx + d.n - 1;
				stride    = d.n;
				break;
            default:
                break;
		}
		coeff = 1.0 / (h * h);
	}
	int row(int i) { return start + stride * i; }
	int size(int i) { return 1; }
	double *coeffs(int i) { return &coeff; }
	int *cols(int i)
	{
		col = nbr_start + stride * i;
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
