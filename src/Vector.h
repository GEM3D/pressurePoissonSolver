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

#ifndef VECTOR_H
#define VECTOR_H
#include <SchurDomain.h>
#include <algorithm>
#include <numeric>
class LocalDataManager
{
	public:
	virtual ~LocalDataManager(){};
};
template <size_t D> class LocalData
{
	private:
	double *                          data;
	std::array<int, D>                strides;
	std::array<int, D>                lengths;
	std::shared_ptr<LocalDataManager> ldm;

	LocalData<D - 1> getSliceOnSidePriv(Side<D> s) const;

	public:
	LocalData(double *data, const std::array<int, D> &strides, const std::array<int, D> &lengths,

	          std::shared_ptr<LocalDataManager> ldm)
	{
		this->data    = data;
		this->strides = strides;
		this->lengths = lengths;
		this->ldm     = ldm;
	}
	inline double &operator[](const std::array<int, D> &coord)
	{
		int idx = 0;
		for (size_t i = 0; i < D; i++) {
			idx += strides[i] * coord[i];
		}
		return data[idx];
	}
	inline const double &operator[](const std::array<int, D> &coord) const
	{
		int idx = 0;
		for (size_t i = 0; i < D; i++) {
			idx += strides[i] * coord[i];
		}
		return data[idx];
	}
	LocalData<D - 1> getSliceOnSide(Side<D> s)
	{
		return getSliceOnSidePriv(s);
	}
	const LocalData<D - 1> getSliceOnSide(Side<D> s) const
	{
		return getSliceOnSidePriv(s);
	}
};
template  <size_t D> inline LocalData<D - 1> LocalData<D>::getSliceOnSidePriv(Side<D> s) const
{
	size_t             axis = s.toInt() / 2;
	std::array<int, D-1> new_strides;
	for (size_t i = 0; i < axis; i++) {
		new_strides[i] = strides[i];
	}
	for (size_t i = axis; i < D-1; i++) {
		new_strides[i] = strides[i + 1];
	}
	std::array<int, D-1> new_lengths;
	for (size_t i = 0; i < axis; i++) {
		new_lengths[i] = lengths[i];
	}
	for (size_t i = axis; i < D-1; i++) {
		new_lengths[i] = lengths[i + 1];
	}
	if (s.isLowerOnAxis()) {
		return LocalData<D - 1>(data, new_strides, new_lengths, ldm);
	} else {
		double *new_data = data + (lengths[axis] - 1) * strides[axis];
		return LocalData<D - 1>(new_data, new_strides, new_lengths, ldm);
	}
}
template <size_t D> class Vector
{
	public:
	virtual ~Vector(){};
	virtual LocalData<D>       getLocalData(int local_patch_id)       = 0;
	virtual const LocalData<D> getLocalData(int local_patch_id) const = 0;
};
#endif
