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
	int n;
	virtual ~StencilHelper() {}
	virtual int     row(int i)    = 0;
	virtual int     size(int i)   = 0;
	virtual double *coeffs(int i) = 0;
	virtual int *   cols(int i)   = 0;
};
class DirichletSH : public StencilHelper
{
	private:
	double coeff;
	int    col;
	int    start;
	int    stride;

	public:
	DirichletSH(PatchInfo<2> &d, Side<2> s)
	{
		double h   = 0;
		int    idx = d.global_index * d.ns[0] * d.ns[1];
		switch (s.toInt()) {
			case Side<2>::west:
				h      = d.spacings[1];
				start  = idx;
				stride = d.ns[0];
				n      = d.ns[1];
				break;
			case Side<2>::east:
				h      = d.spacings[1];
				start  = idx + (d.ns[0] - 1);
				stride = d.ns[0];
				n      = d.ns[1];
				break;
			case Side<2>::south:
				h      = d.spacings[0];
				start  = idx;
				stride = 1;
				n      = d.ns[0];
				break;
			case Side<2>::north:
				h      = d.spacings[0];
				start  = idx + (d.ns[1] - 1) * d.ns[0];
				stride = 1;
				n      = d.ns[0];
				break;
		}
		coeff = -1.0 / (h * h);
	}
	int row(int i)
	{
		return start + stride * i;
	}
	int size(int i)
	{
		return 1;
	}
	double *coeffs(int i)
	{
		return &coeff;
	}
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
	NeumannSH(PatchInfo<2> &d, Side<2> s)
	{
		double h   = 0;
		int    idx = d.global_index * d.ns[0] * d.ns[1];
		switch (s.toInt()) {
			case Side<2>::west:
				h      = d.spacings[1];
				start  = idx;
				stride = d.ns[0];
				n      = d.ns[1];
				break;
			case Side<2>::east:
				h      = d.spacings[1];
				start  = idx + (d.ns[0] - 1);
				stride = d.ns[0];
				n      = d.ns[1];
				break;
			case Side<2>::south:
				h      = d.spacings[0];
				start  = idx;
				stride = 1;
				n      = d.ns[0];
				break;
			case Side<2>::north:
				h      = d.spacings[0];
				start  = idx + (d.ns[1] - 1) * d.ns[0];
				stride = 1;
				n      = d.ns[0];
				break;
		}
		coeff = 1.0 / (h * h);
	}
	int row(int i)
	{
		return start + stride * i;
	}
	int size(int i)
	{
		return 1;
	}
	double *coeffs(int i)
	{
		return &coeff;
	}
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
	NormalSH(PatchInfo<2> &d, Side<2> s)
	{
		double h       = 0;
		int    idx     = d.global_index * d.ns[0] * d.ns[1];
		int    nbr_idx = d.getNormalNbrInfo(s).global_index * d.ns[0] * d.ns[1];
		switch (s.toInt()) {
			case Side<2>::west:
				h         = d.spacings[1];
				start     = idx;
				nbr_start = nbr_idx + d.ns[0] - 1;
				stride    = d.ns[0];
				n         = d.ns[1];
				break;
			case Side<2>::east:
				h         = d.spacings[1];
				start     = idx + (d.ns[0] - 1);
				nbr_start = nbr_idx;
				stride    = d.ns[0];
				n         = d.ns[1];
				break;
			case Side<2>::south:
				h         = d.spacings[0];
				start     = idx;
				nbr_start = nbr_idx + (d.ns[1] - 1) * d.ns[0];
				stride    = 1;
				n         = d.ns[0];
				break;
			case Side<2>::north:
				h         = d.spacings[0];
				start     = idx + (d.ns[1] - 1) * d.ns[0];
				nbr_start = nbr_idx;
				stride    = 1;
				n         = d.ns[0];
				break;
		}
		coeff = 1.0 / (h * h);
	}
	int row(int i)
	{
		return start + stride * i;
	}
	int size(int i)
	{
		return 1;
	}
	double *coeffs(int i)
	{
		return &coeff;
	}
	int *cols(int i)
	{
		col = nbr_start + stride * i;
		return &col;
	}
};
class CoarseSH : public StencilHelper
{
	private:
	std::valarray<double> mid_coeffs = {{-1.0 / 30, -1.0 / 30, 1.0 / 3, 1.0 / 3, 1.0 / 5, 1.0 / 5}};
	std::valarray<double> end_coeffs
	= {{-1.0 / 30, 1.0 / 15, -1.0 / 10, 1.0 / 3, 1.0 / 3, 1.0 / 5, 1.0 / 5}};
	int colz[7];
	int start;
	int bnbr_start;
	int enbr_start;
	int bnbr_start_in;
	int enbr_start_in;
	int stride;

	public:
	CoarseSH(PatchInfo<2> &d, Side<2> s)
	{
		double h             = 0;
		int    idx           = d.global_index * d.ns[0] * d.ns[1];
		int    nbr_idx_left  = d.getFineNbrInfo(s).global_indexes[0] * d.ns[0] * d.ns[1];
		int    nbr_idx_right = d.getFineNbrInfo(s).global_indexes[1] * d.ns[0] * d.ns[1];
		switch (s.toInt()) {
			case Side<2>::west:
				h             = d.spacings[1];
				start         = idx;
				bnbr_start    = nbr_idx_left + d.ns[0] - 1;
				bnbr_start_in = bnbr_start - 1;
				enbr_start    = nbr_idx_right + d.ns[0] - 1;
				enbr_start_in = enbr_start - 1;
				stride        = d.ns[0];
				n             = d.ns[1];
				break;
			case Side<2>::east:
				h             = d.spacings[1];
				start         = idx + (d.ns[0] - 1);
				bnbr_start    = nbr_idx_left;
				bnbr_start_in = bnbr_start + 1;
				enbr_start    = nbr_idx_right;
				enbr_start_in = enbr_start + 1;
				stride        = d.ns[0];
				n             = d.ns[1];
				break;
			case Side<2>::south:
				h             = d.spacings[0];
				start         = idx;
				bnbr_start    = nbr_idx_left + (d.ns[1] - 1) * d.ns[0];
				bnbr_start_in = bnbr_start - d.ns[0];
				enbr_start    = nbr_idx_right + (d.ns[1] - 1) * d.ns[0];
				enbr_start_in = enbr_start - d.ns[0];
				stride        = 1;
				n             = d.ns[0];
				break;
			case Side<2>::north:
				h             = d.spacings[0];
				start         = idx + (d.ns[1] - 1) * d.ns[0];
				bnbr_start    = nbr_idx_left;
				bnbr_start_in = bnbr_start + d.ns[0];
				enbr_start    = nbr_idx_right;
				enbr_start_in = enbr_start + d.ns[0];
				stride        = 1;
				n             = d.ns[0];
				break;
		}
		mid_coeffs /= h * h;
		end_coeffs /= h * h;
	}
	int row(int i)
	{
		return start + stride * i;
	}
	int size(int i)
	{
		if (i == 0 || i == n - 1) {
			return end_coeffs.size();
		} else {
			return mid_coeffs.size();
		}
	}
	double *coeffs(int i)
	{
		if (i == 0 || i == n - 1) {
			return &end_coeffs[0];
		} else {
			return &mid_coeffs[0];
		}
	}
	int *cols(int i)
	{
		if (i == 0) {
			colz[0] = start + stride * 2;
			colz[1] = start + stride * 1;
			colz[2] = start + stride * 0;
			colz[3] = bnbr_start + stride * 0;
			colz[4] = bnbr_start + stride * 1;
			colz[5] = bnbr_start_in + stride * 0;
			colz[6] = bnbr_start_in + stride * 1;
		} else if (i == n - 1) {
			colz[0] = start + stride * (n - 3);
			colz[1] = start + stride * (n - 2);
			colz[2] = start + stride * (n - 1);
			colz[3] = enbr_start + stride * (n - 2);
			colz[4] = enbr_start + stride * (n - 1);
			colz[5] = enbr_start_in + stride * (n - 2);
			colz[6] = enbr_start_in + stride * (n - 1);
		} else {
			colz[0] = start + stride * (i - 1);
			colz[1] = start + stride * (i + 1);
			if (i < n / 2) {
				colz[2] = bnbr_start + stride * 2 * i;
				colz[3] = bnbr_start + stride * (2 * i + 1);
				colz[4] = bnbr_start_in + stride * 2 * i;
				colz[5] = bnbr_start_in + stride * (2 * i + 1);
			} else {
				int j   = i - n / 2;
				colz[2] = enbr_start + stride * 2 * j;
				colz[3] = enbr_start + stride * (2 * j + 1);
				colz[4] = enbr_start_in + stride * 2 * j;
				colz[5] = enbr_start_in + stride * (2 * j + 1);
			}
		}
		return colz;
	}
};
class FineSH : public StencilHelper
{
	private:
	std::valarray<double> end_coeffs = {{1.0 / 12, -3.0 / 10, 3.0 / 4, 2.0 / 3, -1.0 / 5}};
	std::valarray<double> pen_coeffs = {{-1.0 / 20, 7.0 / 30, 7.0 / 20, 2.0 / 3, -1.0 / 5}};
	std::valarray<double> mid_coeffs = {{1.0 / 12, 1.0 / 2, -1.0 / 20, 2.0 / 3, -1.0 / 5}};
	int                   colz[5];
	int                   start;
	int                   start_in;
	int                   nbr_start;
	int                   stride;
	bool                  end;

	public:
	FineSH(PatchInfo<2> &d, Side<2> s)
	{
		double h       = 0;
		int    idx     = d.global_index * d.ns[0] * d.ns[1];
		int    nbr_idx = d.getCoarseNbrInfo(s).global_index * d.ns[0] * d.ns[1];
		end            = d.getCoarseNbrInfo(s).orth_on_coarse.toInt() == 1;
		switch (s.toInt()) {
			case Side<2>::west:
				h         = d.spacings[1];
				start     = idx;
				start_in  = idx + 1;
				nbr_start = nbr_idx + d.ns[0] - 1;
				stride    = d.ns[0];
				n         = d.ns[1];
				break;
			case Side<2>::east:
				h         = d.spacings[1];
				start     = idx + (d.ns[0] - 1);
				start_in  = idx + (d.ns[0] - 2);
				nbr_start = nbr_idx;
				stride    = d.ns[0];
				n         = d.ns[1];
				break;
			case Side<2>::south:
				h         = d.spacings[0];
				start     = idx;
				start_in  = idx + d.ns[0];
				nbr_start = nbr_idx + (d.ns[1] - 1) * d.ns[0];
				stride    = 1;
				n         = d.ns[0];
				break;
			case Side<2>::north:
				h         = d.spacings[0];
				start     = idx + (d.ns[1] - 1) * d.ns[0];
				start_in  = idx + (d.ns[1] - 2) * d.ns[0];
				nbr_start = nbr_idx;
				stride    = 1;
				n         = d.ns[0];
				break;
		}
		end_coeffs /= h * h;
		pen_coeffs /= h * h;
		mid_coeffs /= h * h;
	}
	int row(int i)
	{
		return start + stride * i;
	}
	int size(int i)
	{
		return 5;
	}
	double *coeffs(int i)
	{
		double *retval;
		if (end) {
			if (i == n - 1) {
				retval = &end_coeffs[0];
			} else if (i == n - 2) {
				retval = &pen_coeffs[0];
			} else {
				retval = &mid_coeffs[0];
			}
		} else {
			if (i == 0) {
				retval = &end_coeffs[0];
			} else if (i == 1) {
				retval = &pen_coeffs[0];
			} else {
				retval = &mid_coeffs[0];
			}
		}
		return retval;
	}
	int *cols(int i)
	{
		if (end) {
			if (i == n - 1 || i == n - 2) {
				colz[0] = nbr_start + stride * (n - 3);
				colz[1] = nbr_start + stride * (n - 2);
				colz[2] = nbr_start + stride * (n - 1);
			} else {
				if (i % 2 == 0) {
					colz[0] = nbr_start + stride * (n / 2 + i / 2 - 1);
					colz[1] = nbr_start + stride * (n / 2 + i / 2);
					colz[2] = nbr_start + stride * (n / 2 + i / 2 + 1);
				} else {
					colz[0] = nbr_start + stride * (n / 2 + i / 2 + 1);
					colz[1] = nbr_start + stride * (n / 2 + i / 2);
					colz[2] = nbr_start + stride * (n / 2 + i / 2 - 1);
				}
			}
		} else {
			if (i == 0 || i == 1) {
				colz[0] = nbr_start + stride * 2;
				colz[1] = nbr_start + stride * 1;
				colz[2] = nbr_start + stride * 0;
			} else {
				colz[1] = nbr_start + stride * (i / 2);
				if (i % 2 == 0) {
					colz[0] = nbr_start + stride * (i / 2 - 1);
					colz[2] = nbr_start + stride * (i / 2 + 1);
				} else {
					colz[0] = nbr_start + stride * (i / 2 + 1);
					colz[2] = nbr_start + stride * (i / 2 - 1);
				}
			}
		}
		colz[3] = start + stride * i;
		colz[4] = start_in + stride * i;
		return colz;
	}
};

StencilHelper *getStencilHelper(PatchInfo<2> &d, Side<2> s)
{
	StencilHelper *retval = nullptr;
	if (d.hasNbr(s)) {
		switch (d.getNbrType(s)) {
			case NbrType::Normal:
				retval = new NormalSH(d, s);
				break;
			case NbrType::Coarse:
				retval = new FineSH(d, s);
				break;
			case NbrType::Fine:
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
