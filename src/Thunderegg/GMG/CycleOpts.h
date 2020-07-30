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

#ifndef GMGCYCLEOPTS_H
#define GMGCYCLEOPTS_H
#include <string>
namespace GMG
{
/* case insensitive string
struct ci_char_traits : public char_traits<char> {
    static bool eq(char c1, char c2) { return toupper(c1) == toupper(c2); }
    static bool ne(char c1, char c2) { return toupper(c1) != toupper(c2); }
    static bool lt(char c1, char c2) { return toupper(c1) <  toupper(c2); }
    static int compare(const char* s1, const char* s2, size_t n) {
        while( n-- != 0 ) {
            if( toupper(*s1) < toupper(*s2) ) return -1;
            if( toupper(*s1) > toupper(*s2) ) return 1;
            ++s1; ++s2;
        }
        return 0;
    }
    static const char* find(const char* s, int n, char a) {
        while( n-- > 0 && toupper(*s) != toupper(a) ) {
            ++s;
        }
        return s;
    }
};

typedef std::basic_string<char, ci_char_traits> ci_string;
*/

struct CycleOpts {
	/**
	 * @brief The max number of levels in GMG cycle. 0 means no limit.
	 */
	int max_levels = 0;
	/**
	 * @brief Lowest level is guaranteed to have at least this number of patches per processor.
	 */
	double patches_per_proc = 0;
	/**
	 * @brief Number of sweeps on down cycle
	 */
	int pre_sweeps = 1;
	/**
	 * @brief Number of sweeps on up cycle
	 */
	int post_sweeps = 1;
	/**
	 * @brief Number of sweeps inbetween up and down
	 */
	int mid_sweeps = 1;
	/**
	 * @brief Number of sweeps on coarse level
	 */
	int coarse_sweeps = 1;
	/**
	 * @brief Cycle type
	 */
	std::string cycle_type = "V";
};
} // namespace GMG
#endif
