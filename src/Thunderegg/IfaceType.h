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

#ifndef IFACETYPE_H
#define IFACETYPE_H
#include <Thunderegg/Side.h>
/*
enum class IfaceType {
    normal,
    coarse_to_coarse,
    fine_to_coarse_0,
    fine_to_coarse_1,
    fine_to_coarse_2,
    fine_to_coarse_3,
    fine_to_fine_0,
    fine_to_fine_1,
    fine_to_fine_2,
    fine_to_fine_3,
    coarse_to_fine_0,
    coarse_to_fine_1,
    coarse_to_fine_2,
    coarse_to_fine_3
};
inline IfaceType operator+(const IfaceType &a, const int &b)
{
    return static_cast<IfaceType>(static_cast<int>(a) + b);
}
*/
template <size_t D> class IfaceType
{
	private:
	/**
	 * @brief the value of the enum
	 */
	char           val     = -1;
	Orthant<D - 1> orthant = -1;

	public:
	// enum definitions
	static constexpr char normal           = 0b000;
	static constexpr char coarse_to_coarse = 0b001;
	static constexpr char fine_to_coarse   = 0b010;
	static constexpr char fine_to_fine     = 0b011;
	static constexpr char coarse_to_fine   = 0b100;
	IfaceType()                            = default;
	IfaceType(const char val)
	{
		this->val = val;
	}
	IfaceType(const char val, const Orthant<D - 1> orthant)
	{
		this->val     = val;
		this->orthant = orthant;
	}
	void set(int blah) {}
	char toInt()
	{
		return val;
	}
	Orthant<D - 1> getOrthant()
	{
		return orthant;
	}
	void setOrthant(const Orthant<D - 1> orthant)
	{
		this->orthant = orthant;
	}
	bool operator<(const IfaceType &b) const
	{
		int orth   = orthant.toInt();
		int b_orth = b.orthant.toInt();
		return std::tie(val, orth) < std::tie(b.val, b_orth);
	}
};

#endif
