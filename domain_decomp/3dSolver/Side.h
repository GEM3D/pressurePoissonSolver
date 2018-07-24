#ifndef SIDEENUM_H
#define SIDEENUM_H
#include <array>
#include <iostream>
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
class Side
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
	static std::array<Side, 6> getValues();
	/**
	 * @brief Get an array of all pairs of sides that touch.
	 *
	 * @return The array.
	 */
	static std::array<std::array<Side, 2>, 12> getPairValues();
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
inline std::array<Side, 6> Side::getValues()
{
	return std::array<Side, 6>(
	{{Side::west, Side::east, Side::south, Side::north, Side::bottom, Side::top}});
}
inline std::array<std::array<Side, 2>, 12> Side::getPairValues()
{
	return std::array<std::array<Side, 2>, 12>({{{{Side::west, Side::south}},
	                                             {{Side::west, Side::north}},
	                                             {{Side::west, Side::bottom}},
	                                             {{Side::west, Side::top}},
	                                             {{Side::east, Side::south}},
	                                             {{Side::east, Side::north}},
	                                             {{Side::east, Side::bottom}},
	                                             {{Side::east, Side::top}},
	                                             {{Side::south, Side::bottom}},
	                                             {{Side::south, Side::top}},
	                                             {{Side::north, Side::bottom}},
	                                             {{Side::north, Side::top}}}});
}
inline Side Side::opposite() const
{
	Side retval = *this;
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
class Octant
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

	/**
	 * @brief Default constructor that initializes the value to -1.
	 */
	Octant() = default;
	/**
	 * @brief Create new Octant with given value.
	 *
	 * @param val the value
	 */
	Octant(const int val)
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
	Octant getInteriorNbrOnSide(Side s) const;
	/**
	 * @brief Return the octant that neighbors this octant on a particular side.
	 *
	 * @param s the side of the octant that you want the neighbor of.
	 *
	 * @return  The octant that neighbors on that side.
	 */
	Octant getExteriorNbrOnSide(Side s) const;
	/**
	 * @brief Get the sides of the octant that are on the interior of the cube.
	 *
	 * @return The sides of the octant that are on the interior of the cube.
	 */
	inline std::array<Side, 3> getInteriorSides() const
	{
		std::array<Side, 3> retval;
		retval[0] = val & 0b001 ? Side::west : Side::east;
		retval[1] = val & 0b010 ? Side::south : Side::north;
		retval[2] = val & 0b100 ? Side::bottom : Side::top;
		return retval;
	}
	/**
	 * @brief Get the sides of the octant that are on the exterior of the cube.
	 *
	 * @return The sides of the octant that are on the exterior of the cube.
	 */
	inline std::array<Side, 3> getExteriorSides() const
	{
		std::array<Side, 3> retval;
		retval[0] = val & 0b001 ? Side::east : Side::west;
		retval[1] = val & 0b010 ? Side::north : Side::south;
		retval[2] = val & 0b100 ? Side::top : Side::bottom;
		return retval;
	}
	/**
	 * @brief Return whether or not the octant lies on a particular side of a cube.
	 *
	 * @param s the side of the cube.j
	 *
	 * @return Whether or not it lies on that side.
	 */
	inline bool isOnSide(Side s) const
	{
		int  idx        = s.toInt() / 2;
		int  remainder  = s.toInt() % 2;
		bool is_bit_set = val & (0x1 << idx);
		return is_bit_set == remainder;
	}
	/**
	 * @brief Get an array of all Octant values, in increasing order.
	 *
	 * @return The array.
	 */
	static std::array<Octant, 8> getValues();
	/**
	 * @brief Get an array of all Octant values that lie on a particular side of the cube.
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
	static std::array<Octant, 4> getValuesOnSide(Side s);
	/**
	 * @brief Equals operator.
	 *
	 * @param other The other octant.
	 *
	 * @return Whether or not the value of this octant equals the value other octant.
	 */
	bool operator==(const Octant &other) const
	{
		return val == other.val;
	}
};

// function definitions
inline Octant Octant::getInteriorNbrOnSide(Side s) const
{
	Octant retval = *this;
	// flip the bit for that side
	retval.val ^= (0x1 << (s.toInt() / 2));
	return retval;
}
inline Octant Octant::getExteriorNbrOnSide(Side s) const
{
	Octant retval = *this;
	// flip the bit for that side
	retval.val ^= (0x1 << (s.toInt() / 2));
	return retval;
}

inline std::array<Octant, 4> Octant::getValuesOnSide(Side s)
{
	std::array<Octant, 4> retval;
	switch (s.toInt()) {
		case Side::west:
			retval = {{Octant::bsw, Octant::bnw, Octant::tsw, Octant::tnw}};
			break;
		case Side::east:
			retval = {{Octant::bse, Octant::bne, Octant::tse, Octant::tne}};
			break;
		case Side::south:
			retval = {{Octant::bsw, Octant::bse, Octant::tsw, Octant::tse}};
			break;
		case Side::north:
			retval = {{Octant::bnw, Octant::bne, Octant::tnw, Octant::tne}};
			break;
		case Side::bottom:
			retval = {{Octant::bsw, Octant::bse, Octant::bnw, Octant::bne}};
			break;
		case Side::top:
			retval = {{Octant::tsw, Octant::tse, Octant::tnw, Octant::tne}};
			break;
	}
	return retval;
}
inline std::array<Octant, 8> Octant::getValues()
{
	return std::array<Octant, 8>({{Octant::bsw, Octant::bse, Octant::bnw, Octant::bne, Octant::tsw,
	                               Octant::tse, Octant::tnw, Octant::tne}});
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
inline std::ostream &operator<<(std::ostream &os, const Side &s)
{
	switch (s.toInt()) {
		case Side::west:
			os << "Side::west";
			break;
		case Side::east:
			os << "Side::east";
			break;
		case Side::south:
			os << "Side::south";
			break;
		case Side::north:
			os << "Side::north";
			break;
		case Side::bottom:
			os << "Side::bottom";
			break;
		case Side::top:
			os << "Side::top";
			break;
	}
	return os;
}
#endif
