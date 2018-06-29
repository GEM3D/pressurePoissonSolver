#ifndef STENCILHELPER_H
#define STENCILHELPER_H
#include "Domain.h"
#include <valarray>
class StencilHelper
{
	public:
	virtual ~StencilHelper() {}
	virtual int     row(int xi, int yi)    = 0;
	virtual int     size(int xi, int yi)   = 0;
	virtual double *coeffs(int xi, int yi) = 0;
	virtual int *   cols(int xi, int yi)   = 0;
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
		switch (s.toInt()) {
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
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return 1;
	}
	double *coeffs(int xi, int yi)
	{
		return &coeff;
	}
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
		switch (s.toInt()) {
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
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return 1;
	}
	double *coeffs(int xi, int yi)
	{
		return &coeff;
	}
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
		NormalNbrInfo &nbr_info = d.getNormalNbrInfo(s);
		double         h        = 0;
		int            idx      = d.id_global * d.n * d.n * d.n;
		int            nbr_idx  = nbr_info.global_index * d.n * d.n * d.n;
		switch (s.toInt()) {
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
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return 1;
	}
	double *coeffs(int xi, int yi)
	{
		return &coeff;
	}
	int *cols(int xi, int yi)
	{
		col = nbr_start + stridex * xi + stridey * yi;
		return &col;
	}
};
class CoarseSH : public StencilHelper
{
	private:
	std::valarray<double> coeff = {{5.0 / 6, -1.0 / 6, -1.0 / 6, -1.0 / 6, 4.0 / 6}};
	int                   colz[5];
	int                   quad;
	int                   start;
	int                   nbr_start;
	int                   stridex;
	int                   stridey;
	int                   n;

	public:
	CoarseSH(Domain &d, Side s)
	{
		CoarseNbrInfo &nbr_info = d.getCoarseNbrInfo(s);
		double         h        = 0;
		n                       = d.n;
		int idx                 = d.id_global * d.n * d.n * d.n;
		int nbr_idx             = nbr_info.global_index * d.n * d.n * d.n;
		quad                    = nbr_info.quad_on_coarse;
		switch (s.toInt()) {
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
		coeff /= (h * h);
	}
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return coeff.size();
	}
	double *coeffs(int xi, int yi)
	{
		return &coeff[0];
	}
	int *cols(int xi, int yi)
	{
		switch (quad) {
			case 0:
				colz[4] = nbr_start + stridex * (xi / 2) + stridey * (yi / 2);
				break;
			case 1:
				colz[4] = nbr_start + stridex * ((xi + n) / 2) + stridey * (yi / 2);
				break;
			case 2:
				colz[4] = nbr_start + stridex * (xi / 2) + stridey * ((yi + n) / 2);
				break;
			case 3:
				colz[4] = nbr_start + stridex * ((xi + n) / 2) + stridey * ((yi + n) / 2);
				break;
			default:
				break;
		}
		int nxi;
		if (xi % 2 == 0) {
			nxi = xi + 1;
		} else {
			nxi = xi - 1;
		}
		int nyi = 0;
		if (yi % 2 == 0) {
			nyi = yi + 1;
		} else {
			nyi = yi - 1;
		}
		colz[0] = start + stridex * xi + stridey * yi;
		colz[1] = start + stridex * nxi + stridey * yi;
		colz[2] = start + stridex * xi + stridey * nyi;
		colz[3] = start + stridex * nxi + stridey * nyi;
		return colz;
	}
};
class FineSH : public StencilHelper
{
	private:
	std::valarray<double> coeff = {{-1.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 3, 1.0 / 3}};
	int                   colz[5];
	int                   start;
	int                   nbr_start[4];
	int                   stridex;
	int                   stridey;
	int                   n;

	public:
	FineSH(Domain &d, Side s)
	{
		FineNbrInfo &nbr_info = d.getFineNbrInfo(s);
		double       h        = 0;
		n                     = d.n;
		int idx               = d.id_global * d.n * d.n * d.n;
		int nbr_idx[4];
		for (int i = 0; i < 4; i++) {
			nbr_idx[i] = nbr_info.global_indexes[i] * d.n * d.n * d.n;
		}
		switch (s.toInt()) {
			case Side::west:
				h     = d.x_length / d.n;
				start = idx;
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i] + (d.n - 1);
				}
				stridex = d.n;
				stridey = d.n * d.n;
				break;
			case Side::east:
				h     = d.x_length / d.n;
				start = idx + (d.n - 1);
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i];
				}
				stridex = d.n;
				stridey = d.n * d.n;
				break;
			case Side::south:
				h     = d.y_length / d.n;
				start = idx;
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i] + (d.n - 1) * d.n;
				}
				stridex = 1;
				stridey = d.n * d.n;
				break;
			case Side::north:
				h     = d.y_length / d.n;
				start = idx + (d.n - 1) * d.n;
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i];
				}
				stridex = 1;
				stridey = d.n * d.n;
				break;
			case Side::bottom:
				h     = d.z_length / d.n;
				start = idx;
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i] + (d.n - 1) * d.n * d.n;
				}
				stridex = 1;
				stridey = d.n;
				break;
			case Side::top:
				h     = d.z_length / d.n;
				start = idx + (d.n - 1) * d.n * d.n;
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i];
				}
				stridex = 1;
				stridey = d.n;
				break;
		}
		coeff /= (h * h);
	}
	int row(int xi, int yi)
	{
		return start + stridex * xi + stridey * yi;
	}
	int size(int xi, int yi)
	{
		return coeff.size();
	}
	double *coeffs(int xi, int yi)
	{
		return &coeff[0];
	}
	int *cols(int xi, int yi)
	{
		colz[0]  = start + stridex * xi + stridey * yi;
		int quad = (xi >= n / 2) | ((yi >= n / 2) << 1);
		int nxi  = xi % (n / 2) * 2;
		int nyi  = yi % (n / 2) * 2;
		colz[1]  = nbr_start[quad] + stridex * (nxi) + stridey * (nyi);
		colz[2]  = nbr_start[quad] + stridex * (nxi + 1) + stridey * (nyi);
		colz[3]  = nbr_start[quad] + stridex * (nxi) + stridey * (nyi + 1);
		colz[4]  = nbr_start[quad] + stridex * (nxi + 1) + stridey * (nyi + 1);
		return colz;
	}
};
StencilHelper *getStencilHelper(Domain &d, Side s)
{
	StencilHelper *retval = nullptr;
	if (d.hasNbr(s)) {
		switch (d.getNbrType(s)) {
			case NbrType::Normal:
				retval = new NormalSH(d, s);
				break;
			case NbrType::Fine:
				retval = new FineSH(d, s);
				break;
			case NbrType::Coarse:
				retval = new CoarseSH(d, s);
				break;
		}
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
