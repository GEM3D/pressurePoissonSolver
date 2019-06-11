/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#ifndef STENCILHELPER_H
#define STENCILHELPER_H
#include <Thunderegg/PatchInfo.h>
#include <valarray>
class StencilHelper
{
	public:
	int nx;
	int ny;
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
	DirichletSH(PatchInfo<3> &d, Side<3> s)
	{
		double h   = 0;
		int    idx = d.global_index * d.ns[0] * d.ns[1] * d.ns[2];
		switch (s.toInt()) {
			case Side<3>::west:
				h       = d.spacings[0];
				start   = idx;
				stridex = d.ns[0];
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[1];
				ny      = d.ns[2];
				break;
			case Side<3>::east:
				h       = d.spacings[0];
				start   = idx + (d.ns[0] - 1);
				stridex = d.ns[0];
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[1];
				ny      = d.ns[2];
				break;
			case Side<3>::south:
				h       = d.spacings[1];
				start   = idx;
				stridex = 1;
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[0];
				ny      = d.ns[2];
				break;
			case Side<3>::north:
				h       = d.spacings[1];
				start   = idx + (d.ns[1] - 1) * d.ns[0];
				stridex = 1;
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[0];
				ny      = d.ns[2];
				break;
			case Side<3>::bottom:
				h       = d.spacings[2];
				start   = idx;
				stridex = 1;
				stridey = d.ns[0];
				nx      = d.ns[0];
				ny      = d.ns[1];
				break;
			case Side<3>::top:
				h       = d.spacings[2];
				start   = idx + (d.ns[2] - 1) * d.ns[1] * d.ns[0];
				stridex = 1;
				stridey = d.ns[0];
				nx      = d.ns[0];
				ny      = d.ns[1];
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
	NeumannSH(PatchInfo<3> &d, Side<3> s)
	{
		double h   = 0;
		int    idx = d.global_index * d.ns[0] * d.ns[1] * d.ns[2];
		switch (s.toInt()) {
			case Side<3>::west:
				h       = d.spacings[0];
				start   = idx;
				stridex = d.ns[0];
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[1];
				ny      = d.ns[2];
				break;
			case Side<3>::east:
				h       = d.spacings[0];
				start   = idx + (d.ns[0] - 1);
				stridex = d.ns[0];
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[1];
				ny      = d.ns[2];
				break;
			case Side<3>::south:
				h       = d.spacings[1];
				start   = idx;
				stridex = 1;
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[0];
				ny      = d.ns[2];
				break;
			case Side<3>::north:
				h       = d.spacings[1];
				start   = idx + (d.ns[1] - 1) * d.ns[0];
				stridex = 1;
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[0];
				ny      = d.ns[2];
				break;
			case Side<3>::bottom:
				h       = d.spacings[2];
				start   = idx;
				stridex = 1;
				stridey = d.ns[0];
				nx      = d.ns[0];
				ny      = d.ns[1];
				break;
			case Side<3>::top:
				h       = d.spacings[2];
				start   = idx + (d.ns[2] - 1) * d.ns[1] * d.ns[0];
				stridex = 1;
				stridey = d.ns[0];
				nx      = d.ns[0];
				ny      = d.ns[1];
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
	NormalSH(PatchInfo<3> &d, Side<3> s)
	{
		NormalNbrInfo<3> &nbr_info = d.getNormalNbrInfo(s);
		double            h        = 0;
		int               idx      = d.global_index * d.ns[0] * d.ns[1] * d.ns[2];
		int               nbr_idx  = nbr_info.global_index * d.ns[0] * d.ns[1] * d.ns[2];
		switch (s.toInt()) {
			case Side<3>::west:
				h         = d.spacings[0];
				start     = idx;
				nbr_start = nbr_idx + (d.ns[0] - 1);
				stridex   = d.ns[0];
				stridey   = d.ns[0] * d.ns[1];
				nx        = d.ns[1];
				ny        = d.ns[2];
				break;
			case Side<3>::east:
				h         = d.spacings[0];
				start     = idx + (d.ns[0] - 1);
				nbr_start = nbr_idx;
				stridex   = d.ns[0];
				stridey   = d.ns[0] * d.ns[1];
				nx        = d.ns[1];
				ny        = d.ns[2];
				break;
			case Side<3>::south:
				h         = d.spacings[1];
				start     = idx;
				nbr_start = nbr_idx + (d.ns[1] - 1) * d.ns[0];
				stridex   = 1;
				stridey   = d.ns[0] * d.ns[1];
				nx        = d.ns[0];
				ny        = d.ns[2];
				break;
			case Side<3>::north:
				h         = d.spacings[1];
				start     = idx + (d.ns[1] - 1) * d.ns[0];
				nbr_start = nbr_idx;
				stridex   = 1;
				stridey   = d.ns[0] * d.ns[1];
				nx        = d.ns[0];
				ny        = d.ns[2];
				break;
			case Side<3>::bottom:
				h         = d.spacings[2];
				start     = idx;
				nbr_start = nbr_idx + (d.ns[2] - 1) * d.ns[1] * d.ns[0];
				stridex   = 1;
				stridey   = d.ns[0];
				nx        = d.ns[0];
				ny        = d.ns[1];
				break;
			case Side<3>::top:
				h         = d.spacings[2];
				start     = idx + (d.ns[2] - 1) * d.ns[1] * d.ns[0];
				nbr_start = nbr_idx;
				stridex   = 1;
				stridey   = d.ns[0];
				nx        = d.ns[0];
				ny        = d.ns[1];
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

	public:
	CoarseSH(PatchInfo<3> &d, Side<3> s)
	{
		CoarseNbrInfo<3> &nbr_info = d.getCoarseNbrInfo(s);
		double            h        = 0;
		int               idx      = d.global_index * d.ns[0] * d.ns[1] * d.ns[2];
		int               nbr_idx  = nbr_info.global_index * d.ns[0] * d.ns[1] * d.ns[2];
		quad                       = nbr_info.quad_on_coarse;
		switch (s.toInt()) {
			case Side<3>::west:
				h         = d.spacings[0];
				start     = idx;
				nbr_start = nbr_idx + (d.ns[0] - 1);
				stridex   = d.ns[0];
				stridey   = d.ns[0] * d.ns[1];
				nx        = d.ns[1];
				ny        = d.ns[2];
				break;
			case Side<3>::east:
				h         = d.spacings[0];
				start     = idx + (d.ns[0] - 1);
				nbr_start = nbr_idx;
				stridex   = d.ns[0];
				stridey   = d.ns[0] * d.ns[1];
				nx        = d.ns[1];
				ny        = d.ns[2];
				break;
			case Side<3>::south:
				h         = d.spacings[1];
				start     = idx;
				nbr_start = nbr_idx + (d.ns[1] - 1) * d.ns[0];
				stridex   = 1;
				stridey   = d.ns[0] * d.ns[1];
				nx        = d.ns[0];
				ny        = d.ns[2];
				break;
			case Side<3>::north:
				h         = d.spacings[1];
				start     = idx + (d.ns[1] - 1) * d.ns[0];
				nbr_start = nbr_idx;
				stridex   = 1;
				stridey   = d.ns[0] * d.ns[1];
				nx        = d.ns[0];
				ny        = d.ns[2];
				break;
			case Side<3>::bottom:
				h         = d.spacings[2];
				start     = idx;
				nbr_start = nbr_idx + (d.ns[2] - 1) * d.ns[1] * d.ns[0];
				stridex   = 1;
				stridey   = d.ns[0];
				nx        = d.ns[0];
				ny        = d.ns[1];
				break;
			case Side<3>::top:
				h         = d.spacings[2];
				start     = idx + (d.ns[2] - 1) * d.ns[1] * d.ns[0];
				nbr_start = nbr_idx;
				stridex   = 1;
				stridey   = d.ns[0];
				nx        = d.ns[0];
				ny        = d.ns[1];
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
				colz[4] = nbr_start + stridex * ((xi + nx) / 2) + stridey * (yi / 2);
				break;
			case 2:
				colz[4] = nbr_start + stridex * (xi / 2) + stridey * ((yi + ny) / 2);
				break;
			case 3:
				colz[4] = nbr_start + stridex * ((xi + nx) / 2) + stridey * ((yi + ny) / 2);
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

	public:
	FineSH(PatchInfo<3> &d, Side<3> s)
	{
		FineNbrInfo<3> &nbr_info = d.getFineNbrInfo(s);
		double          h        = 0;
		int             idx      = d.global_index * d.ns[0] * d.ns[1] * d.ns[2];
		int             nbr_idx[4];
		for (int i = 0; i < 4; i++) {
			nbr_idx[i] = nbr_info.global_indexes[i] * d.ns[0] * d.ns[1] * d.ns[2];
		}
		switch (s.toInt()) {
			case Side<3>::west:
				h     = d.spacings[0];
				start = idx;
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i] + (d.ns[0] - 1);
				}
				stridex = d.ns[0];
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[1];
				ny      = d.ns[2];
				break;
			case Side<3>::east:
				h     = d.spacings[0];
				start = idx + (d.ns[0] - 1);
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i];
				}
				stridex = d.ns[0];
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[1];
				ny      = d.ns[2];
				break;
			case Side<3>::south:
				h     = d.spacings[1];
				start = idx;
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i] + (d.ns[1] - 1) * d.ns[0];
				}
				stridex = 1;
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[0];
				ny      = d.ns[2];
				break;
			case Side<3>::north:
				h     = d.spacings[1];
				start = idx + (d.ns[1] - 1) * d.ns[0];
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i];
				}
				stridex = 1;
				stridey = d.ns[0] * d.ns[1];
				nx      = d.ns[0];
				ny      = d.ns[2];
				break;
			case Side<3>::bottom:
				h     = d.spacings[2];
				start = idx;
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i] + (d.ns[2] - 1) * d.ns[1] * d.ns[0];
				}
				stridex = 1;
				stridey = d.ns[0];
				nx      = d.ns[0];
				ny      = d.ns[1];
				break;
			case Side<3>::top:
				h     = d.spacings[2];
				start = idx + (d.ns[2] - 1) * d.ns[1] * d.ns[0];
				for (int i = 0; i < 4; i++) {
					nbr_start[i] = nbr_idx[i];
				}
				stridex = 1;
				stridex = 1;
				stridey = d.ns[0];
				nx      = d.ns[0];
				ny      = d.ns[1];
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
		int quad = (xi >= nx / 2) | ((yi >= ny / 2) << 1);
		int nxi  = xi % (nx / 2) * 2;
		int nyi  = yi % (ny / 2) * 2;
		colz[1]  = nbr_start[quad] + stridex * (nxi) + stridey * (nyi);
		colz[2]  = nbr_start[quad] + stridex * (nxi + 1) + stridey * (nyi);
		colz[3]  = nbr_start[quad] + stridex * (nxi) + stridey * (nyi + 1);
		colz[4]  = nbr_start[quad] + stridex * (nxi + 1) + stridey * (nyi + 1);
		return colz;
	}
};
StencilHelper *getStencilHelper(PatchInfo<3> &d, Side<3> s)
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
