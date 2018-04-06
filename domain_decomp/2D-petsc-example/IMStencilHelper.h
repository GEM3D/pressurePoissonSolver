#ifndef IMSTENCILHELPER_H
#define IMSTENCILHELPER_H
#include "Domain.h"
#include <valarray>
class IMStencilHelper
{
	public:
	virtual ~IMStencilHelper() {}
	virtual int     row(int i)    = 0;
	virtual int     size(int i)   = 0;
	virtual double *coeffs(int i) = 0;
	virtual int *   cols(int i)   = 0;
};
class DirichletSH : public IMStencilHelper
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
		}
		coeff = -1.0 / (h * h);
	}
	int     row(int i) { return start + stride * i; }
	int     size(int i) { return 1; }
	double *coeffs(int i) { return &coeff; }
	int *   cols(int i)
	{
		col = start + stride * i;
		return &col;
	}
};
class NeumannSH : public IMStencilHelper
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
		}
		coeff = 1.0 / (h * h);
	}
	int     row(int i) { return start + stride * i; }
	int     size(int i) { return 1; }
	double *coeffs(int i) { return &coeff; }
	int *   cols(int i)
	{
		col = start + stride * i;
		return &col;
	}
};
class DomainToIfaceSH : public IMStencilHelper
{
	private:
	double coeff[2];
	int    colz[2];
	int    start;
	int    nbr_start;
	int    stride;

	public:
	DomainToIfaceSH(Domain &d, Side s, int iface_start)
	{
		double h   = 0;
		int    idx = d.id_global * d.n * d.n;
		nbr_start  = iface_start + d.globalIndex(s) * d.n;
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
		}
		coeff[0] = -1.0 / (h * h);
		coeff[1] = 2.0 / (h * h);
	}
	int     row(int i) { return start + stride * i; }
	int     size(int i) { return 2; }
	double *coeffs(int i) { return coeff; }
	int *   cols(int i)
	{
		colz[0] = start + stride * i;
		colz[1] = nbr_start + i;
		return colz;
	}
};
class NormalIfaceSH : public IMStencilHelper
{
	private:
	double coeff[3] = {1, -.5, -.5};
	int    colz[3];
	int    start;
	int    left_nbr_start;
	int    right_nbr_start;
	int    stride;

	public:
	NormalIfaceSH(Iface &iface, int iface_start, int n)
	{
		start = iface_start + iface.id_global * n;
		if (iface.y_axis) {
			stride          = n;
			left_nbr_start  = iface.left.id_global * n * n + (n - 1);
			right_nbr_start = iface.right.id_global * n * n;
		} else {
			stride          = 1;
			left_nbr_start  = iface.left.id_global * n * n + n * (n - 1);
			right_nbr_start = iface.right.id_global * n * n;
		}
	}
	int     row(int i) { return start + i; }
	int     size(int i) { return 3; }
	double *coeffs(int i) { return coeff; }
	int *   cols(int i)
	{
		colz[0] = start + i;
		colz[1] = left_nbr_start + stride * i;
		colz[2] = right_nbr_start + stride * i;
		return colz;
	}
};
class CoarseIfaceSH : public IMStencilHelper
{
	private:
	double end_coeff[8]
	= {1, -1.0 / 6, -1.0 / 10, -1.0 / 6, -1.0 / 10, 1.0 / 60, -1.0 / 30, -9.0 / 20};
	double coeff[8] = {1, -1.0 / 6, -1.0 / 10, -1.0 / 6, -1.0 / 10, 1.0 / 60,-.5, 1.0 / 60};
	int    colz[8];
	int    start;
	int    coarse_nbr_start;
	int    left_fine_nbr_start;
	int    left_in_fine_nbr_start;
	int    right_fine_nbr_start;
	int    right_in_fine_nbr_start;
	int    stride;
	int    n;

	public:
	CoarseIfaceSH(Iface &iface, int iface_start, int n)
	{
		this->n = n;
		start   = iface_start + iface.id_global * n;
		if (iface.type == IfaceType::coarse_on_left) {
			if (iface.y_axis) {
				stride                  = n;
				coarse_nbr_start        = iface.left.id_global * n * n + (n - 1);
				left_fine_nbr_start     = iface.extra.id_global * n * n;
				left_in_fine_nbr_start  = iface.extra.id_global * n * n + 1;
				right_fine_nbr_start    = iface.right.id_global * n * n;
				right_in_fine_nbr_start = iface.right.id_global * n * n + 1;
			} else {
				stride                  = 1;
				coarse_nbr_start        = iface.left.id_global * n * n + n * (n - 1);
				left_fine_nbr_start     = iface.right.id_global * n * n;
				left_in_fine_nbr_start  = iface.right.id_global * n * n + n;
				right_fine_nbr_start    = iface.extra.id_global * n * n;
				right_in_fine_nbr_start = iface.extra.id_global * n * n + n;
			}
		} else {
			if (iface.y_axis) {
				stride                  = n;
				left_fine_nbr_start     = iface.left.id_global * n * n + (n - 1);
				left_in_fine_nbr_start  = iface.left.id_global * n * n + (n - 2);
				right_fine_nbr_start    = iface.extra.id_global * n * n + (n - 1);
				right_in_fine_nbr_start = iface.extra.id_global * n * n + (n - 2);
				coarse_nbr_start        = iface.right.id_global * n * n;
			} else {
				stride                  = 1;
				left_fine_nbr_start     = iface.extra.id_global * n * n + n * (n - 1);
				left_in_fine_nbr_start  = iface.extra.id_global * n * n + n * (n - 2);
				right_fine_nbr_start    = iface.left.id_global * n * n + n * (n - 1);
				right_in_fine_nbr_start = iface.left.id_global * n * n + n * (n - 2);
				coarse_nbr_start        = iface.right.id_global * n * n;
			}
		}
	}
	int row(int i) { return start + i; }
	int size(int i)
	{
			return 8;
	}
	double *coeffs(int i)
	{
		if (i != 0 && i != n - 1) {
			return coeff;
		} else {
			return end_coeff;
		}
	}
	int *cols(int i)
	{
		colz[0] = start + i;
		if (i < n / 2) {
			colz[1] = left_fine_nbr_start + (2 * i) * stride;
			colz[2] = left_in_fine_nbr_start + (2 * i) * stride;
			colz[3] = left_fine_nbr_start + (2 * i + 1) * stride;
			colz[4] = left_in_fine_nbr_start + (2 * i + 1) * stride;
		} else {
			colz[1] = right_fine_nbr_start + (2 * (i - n / 2)) * stride;
			colz[2] = right_in_fine_nbr_start + (2 * (i - n / 2)) * stride;
			colz[3] = right_fine_nbr_start + (2 * (i - n / 2) + 1) * stride;
			colz[4] = right_in_fine_nbr_start + (2 * (i - n / 2) + 1) * stride;
		}
		if (i == 0) {
			colz[5] = coarse_nbr_start + (i + 2) * stride;
			colz[6] = coarse_nbr_start + (i + 1) * stride;
			colz[7] = coarse_nbr_start + (i) *stride;
		} else if (i == n - 1) {
			colz[5] = coarse_nbr_start + (i - 2) * stride;
			colz[6] = coarse_nbr_start + (i - 1) * stride;
			colz[7] = coarse_nbr_start + (i) *stride;
		} else {
			colz[5] = coarse_nbr_start + (i - 1) * stride;
			colz[6] = coarse_nbr_start + (i) * stride;
			colz[7] = coarse_nbr_start + (i + 1) * stride;
		}
		return colz;
	}
};
class FineIfaceSH : public IMStencilHelper
{
	private:
	double end_coeff[6] = {1, -1.0 / 24, 3.0 / 20, -3.0 / 8, -5.0 / 6, 1.0 / 10};
	double pen_coeff[6] = {1, 1.0 / 40, -7.0 / 60, -7.0 / 40, -5.0 / 6, 1.0 / 10};
	double coeff[6]     = {1, -1.0 / 24, -1.0 / 4, 1.0 / 40, -5.0 / 6, 1.0 / 10};
	int    colz[6];
	int    start;
	int    coarse_nbr_start;
	int    fine_nbr_start;
	int    in_fine_nbr_start;
	int    stride;
	bool   begin;
	int    n;

	public:
	FineIfaceSH(Iface &iface, int iface_start, int n)
	{
		this->n = n;
		start   = iface_start + iface.id_global * n;
		if (iface.type == IfaceType::refined_on_right_left_of_coarse
		    || iface.type == IfaceType::refined_on_right_right_of_coarse) {
			if (iface.y_axis) {
				stride            = n;
				coarse_nbr_start  = iface.left.id_global * n * n + (n - 1);
				fine_nbr_start    = iface.right.id_global * n * n;
				in_fine_nbr_start = iface.right.id_global * n * n + 1;
			} else {
				stride            = 1;
				coarse_nbr_start  = iface.left.id_global * n * n + n * (n - 1);
				fine_nbr_start    = iface.right.id_global * n * n;
				in_fine_nbr_start = iface.right.id_global * n * n + n;
			}
		} else {
			if (iface.y_axis) {
				stride            = n;
				fine_nbr_start    = iface.left.id_global * n * n + (n - 1);
				in_fine_nbr_start = iface.left.id_global * n * n + (n - 2);
				coarse_nbr_start  = iface.right.id_global * n * n;
			} else {
				stride            = 1;
				fine_nbr_start    = iface.left.id_global * n * n + n * (n - 1);
				in_fine_nbr_start = iface.left.id_global * n * n + n * (n - 2);
				coarse_nbr_start  = iface.right.id_global * n * n;
			}
		}
		if (iface.y_axis) {
			switch (iface.type) {
				case IfaceType::refined_on_right_right_of_coarse:
				case IfaceType::refined_on_left_left_of_coarse:
					begin = true;
					break;
				case IfaceType::refined_on_left_right_of_coarse:
				case IfaceType::refined_on_right_left_of_coarse:
					begin = false;
					break;
				default:
					break;
			}
		} else {
			switch (iface.type) {
				case IfaceType::refined_on_right_right_of_coarse:
				case IfaceType::refined_on_left_left_of_coarse:
					begin = false;
					break;
				case IfaceType::refined_on_left_right_of_coarse:
				case IfaceType::refined_on_right_left_of_coarse:
					begin = true;
					break;
				default:
					break;
			}
		}
	}
	int     row(int i) { return start + i; }
	int     size(int i) { return 6; }
	double *coeffs(int i)
	{
		if (begin) {
			if (i == 0) {
				return end_coeff;
			} else if (i == 1) {
				return pen_coeff;
			} else {
				return coeff;
			}
		} else {
			if (i == n - 1) {
				return end_coeff;
			} else if (i == n - 2) {
				return pen_coeff;
			} else {
				return coeff;
			}
		}
	}
	int *cols(int i)
	{
		colz[0] = start + i;
		if (begin) {
			if (i == 0 || i == 1) {
				colz[1] = coarse_nbr_start + 2 * stride;
				colz[2] = coarse_nbr_start + 1 * stride;
				colz[3] = coarse_nbr_start + 0 * stride;
			} else if (i % 2 == 0) {
				colz[1] = coarse_nbr_start + (i / 2 - 1) * stride;
				colz[2] = coarse_nbr_start + (i / 2) * stride;
				colz[3] = coarse_nbr_start + (i / 2 + 1) * stride;
			} else {
				colz[1] = coarse_nbr_start + (i / 2 + 1) * stride;
				colz[2] = coarse_nbr_start + (i / 2) * stride;
				colz[3] = coarse_nbr_start + (i / 2 - 1) * stride;
			}
		} else {
			if (i == n - 1 || i == n - 2) {
				colz[1] = coarse_nbr_start + (n - 3) * stride;
				colz[2] = coarse_nbr_start + (n - 2) * stride;
				colz[3] = coarse_nbr_start + (n - 1) * stride;
			} else if (i % 2 == 0) {
				colz[1] = coarse_nbr_start + (n / 2 + i / 2 - 1) * stride;
				colz[2] = coarse_nbr_start + (n / 2 + i / 2) * stride;
				colz[3] = coarse_nbr_start + (n / 2 + i / 2 + 1) * stride;
			} else {
				colz[1] = coarse_nbr_start + (n / 2 + i / 2 + 1) * stride;
				colz[2] = coarse_nbr_start + (n / 2 + i / 2) * stride;
				colz[3] = coarse_nbr_start + (n / 2 + i / 2 - 1) * stride;
			}
		}
		colz[4] = fine_nbr_start + i * stride;
		colz[5] = in_fine_nbr_start + i * stride;
		return colz;
	}
};

IMStencilHelper *getIMStencilHelper(Domain &d, Side s, int iface_start)
{
	IMStencilHelper *retval = nullptr;
	if (d.hasNbr(s)) {
		retval = new DomainToIfaceSH(d, s, iface_start);
	} else {
		if (d.isNeumann(s)) {
			retval = new NeumannSH(d, s);
		} else {
			retval = new DirichletSH(d, s);
		}
	}
	return retval;
}
IMStencilHelper *getIfaceStencilHelper(Iface &iface, int iface_start, int n)
{
	IMStencilHelper *retval = nullptr;
	switch (iface.type) {
		case IfaceType::normal:
			retval = new NormalIfaceSH(iface, iface_start, n);
			break;
		case IfaceType::coarse_on_left:
		case IfaceType::coarse_on_right:
			retval = new CoarseIfaceSH(iface, iface_start, n);
			break;
		case IfaceType::refined_on_left_left_of_coarse:
		case IfaceType::refined_on_left_right_of_coarse:
		case IfaceType::refined_on_right_left_of_coarse:
		case IfaceType::refined_on_right_right_of_coarse:
			retval = new FineIfaceSH(iface, iface_start, n);
			break;
	}
	return retval;
}
#endif
