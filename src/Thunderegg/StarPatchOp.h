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

#ifndef STARPATCHOP_H
#define STARPATCHOP_H
#include <Thunderegg/PatchOperator.h>
template <size_t D> class StarPatchOp : public PatchOperator<D>
{
	public:
	void apply(SchurInfo<D> &sinfo, std::shared_ptr<const Vector<D>> u,
	           std::shared_ptr<const Vector<D - 1>> gamma, std::shared_ptr<Vector<D>> f)
	{
		LocalData<D>          f_data = f->getLocalData(sinfo.pinfo->local_index);
		const LocalData<D>    u_data = u->getLocalData(sinfo.pinfo->local_index);
		std::array<double, D> h2     = sinfo.pinfo->spacings;
		for (size_t i = 0; i < D; i++) {
			h2[i] *= h2[i];
		}
		{
			constexpr int axis       = 0;
			Side<D>       lower_side = axis * 2;
			Side<D>       upper_side = axis * 2 + 1;
			if (sinfo.pinfo->hasNbr(lower_side)) {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> bnd
				= gamma->getLocalData(sinfo.getIfaceLocalIndex(lower_side));
				const LocalData<D - 1> mid   = u_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper = u_data.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (2 * bnd[coord] - 3 * mid[coord] + upper[coord]) / h2[axis];
				});
			} else if (sinfo.pinfo->isNeumann(lower_side)) {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u_data.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (-mid[coord] + upper[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u_data.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (-3 * mid[coord] + upper[coord]) / h2[axis];
				});
			}
			// middle
			{
				std::array<int, D> start = f_data.getStart();
				std::array<int, D> end   = f_data.getEnd();
				start[axis] += 1;
				end[axis] -= 1;
				int stride = u_data.getStrides()[axis];
				nested_loop<D>(start, end, [&](std::array<int, D> coord) {
					const double *ptr   = u_data.getPtr(coord);
					double        lower = *(ptr - stride);
					double        mid   = *ptr;
					double        upper = *(ptr + stride);
					f_data[coord]       = (lower - 2 * mid + upper) / h2[axis];
				});
			}
			// east
			if (sinfo.pinfo->hasNbr(upper_side)) {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u_data.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(upper_side);
				const LocalData<D - 1> bnd
				= gamma->getLocalData(sinfo.getIfaceLocalIndex(upper_side));

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (lower[coord] - 3 * mid[coord] + 2 * bnd[coord]) / h2[axis];
				});
			} else if (sinfo.pinfo->isNeumann(upper_side)) {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u_data.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (lower[coord] - mid[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u_data.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] = (lower[coord] - 3 * mid[coord]) / h2[axis];
				});
			}
		}
		for (size_t axis = 1; axis < D; axis++) {
			Side<D> lower_side = axis * 2;
			Side<D> upper_side = axis * 2 + 1;
			if (sinfo.pinfo->hasNbr(lower_side)) {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> bnd
				= gamma->getLocalData(sinfo.getIfaceLocalIndex(lower_side));
				const LocalData<D - 1> mid   = u_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper = u_data.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (2 * bnd[coord] - 3 * mid[coord] + upper[coord]) / h2[axis];
				});
			} else if (sinfo.pinfo->isNeumann(lower_side)) {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u_data.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (-mid[coord] + upper[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(lower_side);
				const LocalData<D - 1> upper   = u_data.getSliceOnSide(lower_side, 1);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (-3 * mid[coord] + upper[coord]) / h2[axis];
				});
			}
			// middle
			{
				std::array<int, D> start = f_data.getStart();
				std::array<int, D> end   = f_data.getEnd();
				start[axis] += 1;
				end[axis] -= 1;
				int stride = u_data.getStrides()[axis];
				nested_loop<D>(start, end, [&](std::array<int, D> coord) {
					const double *ptr   = u_data.getPtr(coord);
					double        lower = *(ptr - stride);
					double        mid   = *ptr;
					double        upper = *(ptr + stride);
					f_data[coord] += (lower - 2 * mid + upper) / h2[axis];
				});
			}
			// east
			if (sinfo.pinfo->hasNbr(upper_side)) {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u_data.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(upper_side);
				const LocalData<D - 1> bnd
				= gamma->getLocalData(sinfo.getIfaceLocalIndex(upper_side));

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (lower[coord] - 3 * mid[coord] + 2 * bnd[coord]) / h2[axis];
				});
			} else if (sinfo.pinfo->isNeumann(upper_side)) {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u_data.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (lower[coord] - mid[coord]) / h2[axis];
				});
			} else {
				LocalData<D - 1>       f_slice = f_data.getSliceOnSide(upper_side);
				const LocalData<D - 1> lower   = u_data.getSliceOnSide(upper_side, 1);
				const LocalData<D - 1> mid     = u_data.getSliceOnSide(upper_side);

				nested_loop<D - 1>(mid.getStart(), mid.getEnd(), [&](std::array<int, D - 1> coord) {
					f_slice[coord] += (lower[coord] - 3 * mid[coord]) / h2[axis];
				});
			}
		}
	}
};
extern template class StarPatchOp<2>;
extern template class StarPatchOp<3>;
#endif
