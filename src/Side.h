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

#ifndef SIDEENUM_H
#define SIDEENUM_H
#include <array>
#include <iostream>
#include <numeric>
#include <vector>
enum class Axis { x, y, z };
/**
 * @brief An enum-style class that represents the sides of a cube.
 *
 * The sides of the cube are named in the following way:
 *
 * Orthogonal to axis | Lower on axis | Higher on axis
 * ------------------ | ------------- | --------------
 *  x-axis            | west          | east
 *  y-axis            | south         | north
 *  z-axis            | bottom        | top
 *
 */
template <size_t D> class Side
{
	private:
	/**
	 * @brief the value of the enum
	 */
	int val = -1;

	public:
	// enum definitions
	static constexpr int west   = 0b000;
	static constexpr int east   = 0b001;
	static constexpr int south  = 0b010;
	static constexpr int north  = 0b011;
	static constexpr int bottom = 0b100;
	static constexpr int top    = 0b101;

	static constexpr int num_sides = 2 * D;

	/**
	 * @brief Default constructor that initializes the value to -1
	 */
	Side() = default;
	/**
	 * @brief Initialize new Side with given value
	 *
	 * @param val the value
	 */
	Side(const int val)
	{
		this->val = val;
	}
	/**
	 * @brief Get the integer value of the side.
	 *
	 * @return The value of the side.
	 */
	int toInt() const
	{
		return val;
	}
	/**
	 * @brief Get an array of all side values, in increasing order.
	 *
	 * @return The array.
	 */
	static std::array<Side<D>, num_sides> getValues();
	/**
	 * @brief Return whether or not the side of the cube is lower on the axis that is orthogonal to
	 * it.
	 *
	 * For example: The x-axis is orthogonal to both the west and east sides. Since west is lower on
	 * the axis, west will return true and east will return false.
	 *
	 * @return Whether or not it is lower on the axis.
	 */
	inline bool isLowerOnAxis() const
	{
		// is least-significant bit set?
		return !(val & 0x1);
	}
	/**
	 * @brief Return the opposite side of the cube.
	 *
	 * For example: the opposite of east is west.
	 *
	 * @return The opposite side.
	 */
	Side opposite() const;
	/**
	 * @brief Compare the enum values.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the value of this side is lower than the other side.
	 */
	bool operator<(const Side &other) const
	{
		return val < other.val;
	}
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other side.
	 *
	 * @return Whether or not the value of this side equals the value other side.
	 */
	bool operator==(const Side &other) const
	{
		return val == other.val;
	}
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other value.
	 *
	 * @return Whether or not the value of this side equals the value of the integer.
	 */
	bool operator==(const int &other) const
	{
		return val == other;
	}
};
template <size_t D> inline std::array<Side<D>, Side<D>::num_sides> Side<D>::getValues()
{
	std::array<Side<D>, Side<D>::num_sides> retval;
	std::iota(retval.begin(), retval.end(), 0);
	return retval;
}
template <size_t D> inline Side<D> Side<D>::opposite() const
{
	Side<D> retval = *this;
	retval.val ^= 0x1;
	return retval;
}
/**
 * @brief An enum-style class that represents the octants of a cube.
 *
 * The octants are named in the following way:
 *
 * bse means Bottom-South-West, which is in the corner of where the bottom, south, and west sides
 * meet.
 */
template <size_t D> class Orthant
{
	private:
	/**
	 * @brief the value of the enum.
	 */
	int val = -1;

	public:
	/**
	 * @brief Bottom-South-West octant of cube.
	 */
	static constexpr int bsw = 0b000;
	/**
	 * @brief Bottom-South-East octant of cube.
	 */
	static constexpr int bse = 0b001;
	/**
	 * @brief Bottom-North-West octant of cube.
	 */
	static constexpr int bnw = 0b010;
	/**
	 * @brief Bottom-North-East octant of cube.
	 */
	static constexpr int bne = 0b011;
	/**
	 * @brief Top-South-West octant of cube.
	 */
	static constexpr int tsw = 0b100;
	/**
	 * @brief Top-South-East octant of cube.
	 */
	static constexpr int tse = 0b101;
	/**
	 * @brief Top-North-West octant of cube.
	 */
	static constexpr int tnw = 0b110;
	/**
	 * @brief Top-North-East octant of cube.
	 */
	static constexpr int tne = 0b111;

	static constexpr int num_orthants = 1 << D;
	/**
	 * @brief Default constructor that initializes the value to -1.
	 */
	Orthant() = default;
	/**
	 * @brief Create new Orthant<D> with given value.
	 *
	 * @param val the value
	 */
	Orthant(const int val)
	{
		this->val = val;
	}
	/**
	 * @brief Get the integer value of the octant.
	 *
	 * @return The integer value.
	 */
	int toInt() const
	{
		return val;
	}
	/**
	 * @brief Return the octant that neighbors this octant on a particular side.
	 *
	 * @param s the side of the octant that you want the neighbor of.
	 *
	 * @return  The octant that neighbors on that side.
	 */
	Orthant<D> getInteriorNbrOnSide(Side<D> s) const;
	/**
	 * @brief Return the octant that neighbors this octant on a particular side.
	 *
	 * @param s the side of the octant that you want the neighbor of.
	 *
	 * @return  The octant that neighbors on that side.
	 */
	Orthant<D> getExteriorNbrOnSide(Side<D> s) const;
	/**
	 * @brief Get the sides of the octant that are on the interior of the cube.
	 *
	 * @return The sides of the octant that are on the interior of the cube.
	 */
	inline std::array<Side<D>, D> getInteriorSides() const
	{
		std::array<Side<D>, D> retval;
		for (size_t i = 0; i < D; i++) {
			int side = 2 * i;
			if (!((1 << i) & val)) { side |= 1; }
			retval[i] = side;
		}
		return retval;
	}
	/**
	 * @brief Get the sides of the octant that are on the exterior of the cube.
	 *
	 * @return The sides of the octant that are on the exterior of the cube.
	 */
	inline std::array<Side<D>, D> getExteriorSides() const
	{
		std::array<Side<D>, D> retval;
		for (size_t i = 0; i < D; i++) {
			int side = 2 * i;
			if ((1 << i) & val) { side |= 1; }
			retval[i] = side;
        }
		return retval;
	}
	/**
	 * @brief Return whether or not the octant lies on a particular side of a cube.
	 *
	 * @param s the side of the cube.j
	 *
	 * @return Whether or not it lies on that side.
	 */
	inline bool isOnSide(Side<D> s) const
	{
		int  idx        = s.toInt() / 2;
		int  remainder  = s.toInt() % 2;
		bool is_bit_set = val & (0x1 << idx);
		return is_bit_set == remainder;
	}
	/**
	 * @brief Get an array of all Orthant<D> values, in increasing order.
	 *
	 * @return The array.
	 */
	static std::array<Orthant, num_orthants> getValues();
	/**
	 * @brief Get an array of all Orthant<D> values that lie on a particular side of the cube.
	 *
	 * When the two axis that the side lies on are arranged in the following way, the octants are
	 * returned in the following order:
	 *
	 *   ^
	 *   |
	 *   |  2  |  3
	 *   |-----+-----
	 *   |  0  |  1
	 *   +----------->
	 *
	 * @return The array.
	 */
	static std::array<Orthant, num_orthants / 2> getValuesOnSide(Side<D> s);
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other octant.
	 *
	 * @return Whether or not the value of this octant equals the value other octant.
	 */
	bool operator==(const Orthant<D> &other) const
	{
		return val == other.val;
	}
};

// function definitions
template <size_t D> inline Orthant<D> Orthant<D>::getInteriorNbrOnSide(Side<D> s) const
{
	Orthant<D> retval = *this;
	// flip the bit for that side
	retval.val ^= (0x1 << (s.toInt() / 2));
	return retval;
}
template <size_t D> inline Orthant<D> Orthant<D>::getExteriorNbrOnSide(Side<D> s) const
{
	Orthant<D> retval = *this;
	// flip the bit for that side
	retval.val ^= (0x1 << (s.toInt() / 2));
	return retval;
}
template <size_t D>
inline std::array<Orthant<D>, Orthant<D>::num_orthants / 2> Orthant<D>::getValuesOnSide(Side<D> s)
{
	size_t bit_to_insert = s.toInt() / 2;
	size_t set_bit       = s.isLowerOnAxis() ? 0 : 1;
	size_t lower_mask    = ~((~0x0) << bit_to_insert);
	size_t upper_mask    = (~0x0) << (bit_to_insert + 1);

	std::array<Orthant<D>, Orthant<D>::num_orthants / 2> retval;
	for (size_t i = 0; i < Orthant<D>::num_orthants / 2; i++) {
		size_t value = (i << 1) & upper_mask;
		value |= i & lower_mask;
		value |= set_bit << bit_to_insert;
		retval[i] = value;
	}
	return retval;
}
template <size_t D> inline std::array<Orthant<D>, Orthant<D>::num_orthants> Orthant<D>::getValues()
{
	std::array<Orthant<D>, Orthant<D>::num_orthants> retval;
	std::iota(retval.begin(), retval.end(), 0);
	return retval;
}
/**
 * @brief ostream operator that prints a string representation of side enum.
 *
 * For example, Side::west will print out "Side::west".
 *
 * @param os the ostream
 * @param s the side to print out.
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Side<2> &s)
{
	switch (s.toInt()) {
		case Side<2>::west:
			os << "Side::west";
			break;
		case Side<2>::east:
			os << "Side::east";
			break;
		case Side<2>::south:
			os << "Side::south";
			break;
		case Side<2>::north:
			os << "Side::north";
			break;
	}
	return os;
}
/**
 * @brief ostream operator that prints a string representation of side enum.
 *
 * For example, Side::west will print out "Side::west".
 *
 * @param os the ostream
 * @param s the side to print out.
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Side<3> &s)
{
	switch (s.toInt()) {
		case Side<3>::west:
			os << "Side::west";
			break;
		case Side<3>::east:
			os << "Side::east";
			break;
		case Side<3>::south:
			os << "Side::south";
			break;
		case Side<3>::north:
			os << "Side::north";
			break;
		case Side<3>::bottom:
			os << "Side::bottom";
			break;
		case Side<3>::top:
			os << "Side::top";
			break;
	}
	return os;
}
/**
 * @brief ostream operator that prints a string representation of orthant enum.
 *
 * For example, Side::west will print out "Orthant::bsw".
 *
 * @param os the ostream
 * @param s the side to print out.
 *
 * @return  the ostream
 */
inline std::ostream &operator<<(std::ostream &os, const Orthant<3> &o)
{
	switch (o.toInt()) {
		case Orthant<3>::bsw:
			os << "Orthant::bsw";
			break;
		case Orthant<3>::bse:
			os << "Orthant::bse";
			break;
		case Orthant<3>::bnw:
			os << "Orthant::bnw";
			break;
		case Orthant<3>::bne:
			os << "Orthant::bne";
			break;
		case Orthant<3>::tsw:
			os << "Orthant::tsw";
			break;
		case Orthant<3>::tse:
			os << "Orthant::tse";
			break;
		case Orthant<3>::tnw:
			os << "Orthant::tnw";
			break;
		case Orthant<3>::tne:
			os << "Orthant::tne";
			break;
	}
	return os;
}
#endif
