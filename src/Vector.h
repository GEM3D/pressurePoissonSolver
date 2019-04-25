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
#include <cmath>
#include <mpi.h>
#include <numeric>
class LocalDataManager
{
	public:
	virtual ~LocalDataManager(){};
};
template <size_t D, size_t Dir, typename T> class NestedLoop
{
	public:
	static void inline nested_loop_loop(std::array<int, D> &coord, std::array<int, D> &start,
	                                    std::array<int, D> &end, T lambda)
	{
		for (coord[Dir] = start[Dir]; coord[Dir] <= end[Dir]; coord[Dir]++) {
			NestedLoop<D, Dir - 1, T>::nested_loop_loop(coord, start, end, lambda);
		}
	}
};

template <size_t D, typename T> class NestedLoop<D, 0, T>
{
	public:
	static void inline nested_loop_loop(std::array<int, D> &coord, std::array<int, D> &start,
	                                    std::array<int, D> &end, T lambda)
	{
		for (coord[0] = start[0]; coord[0] <= end[0]; coord[0]++) {
			lambda(coord);
		}
	}
};
template <size_t D, typename T>
inline void nested_loop(std::array<int, D> start, std::array<int, D> end, T lambda)
{
	std::array<int, D> coord = start;
	NestedLoop<D, D - 1, T>::nested_loop_loop(coord, start, end, lambda);
}
template <size_t D> class LocalData
{
	private:
	double *                          data;
	std::array<int, D>                strides;
	std::array<int, D>                lengths;
	std::array<int, D>                start;
	std::array<int, D>                end;
	std::shared_ptr<LocalDataManager> ldm;

	LocalData<D - 1> getSliceOnSidePriv(Side<D> s) const;

	public:
	LocalData() = default;
	LocalData(double *data, const std::array<int, D> &strides, const std::array<int, D> &lengths,

	          std::shared_ptr<LocalDataManager> ldm)
	{
		this->data    = data;
		this->strides = strides;
		this->lengths = lengths;
		this->ldm     = ldm;
		start.fill(0);
		end = lengths;
		for (size_t i = 0; i < D; i++) {
			end[i]--;
		}
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
	const std::array<int, D> &getLengths() const
	{
		return lengths;
	}
	const std::array<int, D> &getStrides() const
	{
		return strides;
	}
	const std::array<int, D> &getStart() const
	{
		return start;
	}
	const std::array<int, D> &getEnd() const
	{
		return end;
	}
	double *getPtr() const
	{
		return data;
	}
};
template <size_t D> inline LocalData<D - 1> LocalData<D>::getSliceOnSidePriv(Side<D> s) const
{
	size_t                 axis = s.toInt() / 2;
	std::array<int, D - 1> new_strides;
	for (size_t i = 0; i < axis; i++) {
		new_strides[i] = strides[i];
	}
	for (size_t i = axis; i < D - 1; i++) {
		new_strides[i] = strides[i + 1];
	}
	std::array<int, D - 1> new_lengths;
	for (size_t i = 0; i < axis; i++) {
		new_lengths[i] = lengths[i];
	}
	for (size_t i = axis; i < D - 1; i++) {
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
	protected:
	int      num_local_patches;
	MPI_Comm comm = MPI_COMM_WORLD;

	public:
	virtual ~Vector(){};
	virtual LocalData<D>       getLocalData(int local_patch_id)       = 0;
	virtual const LocalData<D> getLocalData(int local_patch_id) const = 0;

	virtual void set(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] = alpha; });
		}
	}
	virtual void scale(double alpha)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] *= alpha; });
		}
	}
	virtual void shift(double delta)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] += delta; });
		}
	}
	virtual void add(std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] += ld_b[coord]; });
		}
	}
	virtual void addScaled(double alpha, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { ld[coord] += ld_b[coord] * alpha; });
		}
	}
	virtual void scaleThenAdd(double alpha, std::shared_ptr<const Vector<D>> b)
	{
		for (int i = 0; i < num_local_patches; i++) {
			LocalData<D>       ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(), [&](std::array<int, D> coord) {
				ld[coord] = alpha * ld[coord] + ld_b[coord];
			});
		}
	}
	virtual double twoNorm() const
	{
		double sum = 0;
		for (int i = 0; i < num_local_patches; i++) {
			const LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { sum += ld[coord] * ld[coord]; });
		}
        MPI_Allreduce(&sum,&sum,1,MPI_DOUBLE,MPI_SUM,comm);
		return sqrt(sum);
	}
	virtual double infNorm() const
	{
		double max = 0;
		for (int i = 0; i < num_local_patches; i++) {
			const LocalData<D> ld = getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { max = fmax(abs(ld[coord]), max); });
		}
        MPI_Allreduce(&max,&max,1,MPI_DOUBLE,MPI_MAX,comm);
		return max;
	}
	virtual double dot(std::shared_ptr<const Vector<D>> b) const
	{
		double retval = 0;
		for (int i = 0; i < num_local_patches; i++) {
			const LocalData<D> ld   = getLocalData(i);
			const LocalData<D> ld_b = b->getLocalData(i);
			nested_loop<D>(ld.getStart(), ld.getEnd(),
			               [&](std::array<int, D> coord) { retval += ld[coord] * ld_b[coord]; });
		}
        MPI_Allreduce(&retval,&retval,1,MPI_DOUBLE,MPI_SUM,comm);
		return retval;
	}
};
extern template class LocalData<1>;
extern template class LocalData<2>;
extern template class LocalData<3>;
extern template class Vector<1>;
extern template class Vector<2>;
extern template class Vector<3>;
#endif
