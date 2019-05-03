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

#ifndef UTILS_H
#define UTILS_H
#include <Thunderegg/SchurDomain.h>
#include <Thunderegg/Vector.h>
#include <array>
#include <cmath>
#include <numeric>
namespace Utils
{
inline int index(const int &n, const int &xi, const int &yi, const int &zi)
{
	return xi + yi * n + zi * n * n;
}
inline int index(const int &n, const int &xi, const int &yi)
{
	return xi + yi * n;
}
template <size_t D> class Slice
{
	private:
	double *           start;
	std::array<int, D> strides;

	public:
	inline Slice() {}
	inline Slice(double *start, std::array<int, D> strides)
	{
		this->start   = start;
		this->strides = strides;
	}
	inline double &operator()(std::array<int, D> coord)
	{
		return start[std::inner_product(strides.begin(), strides.end(), coord.begin(), 0)];
	}
};
template <size_t D> inline Slice<D> getSlice(double *u_view, int n, Side<D + 1> s)
{
	Slice<D>           retval;
	std::array<int, D> strides;
	for (int i = 0; i < s.toInt() / 2; i++) {
		strides[i] = std::pow(n, i);
	}
	for (int i = s.toInt() / 2; i < (int) D; i++) {
		strides[i] = std::pow(n, i + 1);
	}
	if (s.isLowerOnAxis()) {
		retval = Slice<D>(&u_view[0], strides);
	} else {
		retval = Slice<D>(&u_view[(n - 1) * (int) std::pow(n, s.toInt() / 2)], strides);
	}
	return retval;
}
template <size_t D> inline Slice<D> getSlice(double *u_view, int n, Side<D + 1> s, const int *npow)
{
	Slice<D>           retval;
	std::array<int, D> strides;
	for (int i = 0; i < s.toInt() / 2; i++) {
		strides[i] = npow[i];
	}
	for (int i = s.toInt() / 2; i < (int) D; i++) {
		strides[i] = npow[i + 1];
	}
	if (s.isLowerOnAxis()) {
		retval = Slice<D>(&u_view[0], strides);
	} else {
		retval = Slice<D>(&u_view[(n - 1) * npow[s.toInt() / 2]], strides);
	}
	return retval;
}
template <size_t D> inline Slice<D - 1> getSlice(SchurDomain<D> &d, double *u_view, Side<D> s)
{
	int start = d.local_index * std::pow(d.n, D);
	return getSlice<D - 1>(&u_view[start], d.n, s);
}
template <size_t D>
inline Slice<D - 1> getSlice(SchurDomain<D> &d, double *u_view, Side<D> s, const int *npow)
{
	int start = d.local_index * npow[D];
	return getSlice<D - 1>(&u_view[start], d.n, s, npow);
}
} // namespace Utils
#endif
