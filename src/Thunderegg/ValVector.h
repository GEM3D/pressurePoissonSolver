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

#ifndef VALVECTOR_H
#define VALVECTOR_H
#include <Thunderegg/Vector.h>
#include <valarray>
template <size_t D> class ValVector : public Vector<D>

{
	private:
	int                patch_stride;
	std::array<int, D> lengths;
	std::array<int, D> strides;

	public:
	std::valarray<double> vec;
	ValVector() = default;
	ValVector(const std::array<int, D> &lengths, int num_patches = 1)
	{
		int size      = 1;
		this->lengths = lengths;
		for (size_t i = 0; i < D; i++) {
			strides[i] = size;
			size *= lengths[i];
		}
		if (num_patches == 1) {
			patch_stride = 0;
		} else {
			patch_stride = size;
		}
		size *= num_patches;
		vec.resize(size);
	}
	~ValVector() = default;
	LocalData<D> getLocalData(int local_patch_id)
	{
		double *data = &vec[patch_stride * local_patch_id];
		return LocalData<D>(data, strides, lengths, nullptr);
	}
	const LocalData<D> getLocalData(int local_patch_id) const
	{
		double *data = const_cast<double *>(&vec[patch_stride * local_patch_id]);
		return LocalData<D>(data, strides, lengths, nullptr);
	}
};
#endif
