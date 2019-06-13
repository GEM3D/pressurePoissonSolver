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
/**
 * \file
 *
 * @brief Useful type definitions for Thunderegg library
 */
#ifndef THUNDEREGG_TYPEDEFS_H
#define THUNDEREGG_TYPEDEFS_H

#include <array>
#include <functional>

/**
 * @brief This function is used to set the boundary condition types of the solver.
 *
 * This function will only be called on the boundaries of the domain.
 *
 * @param s the side of the patch being inquired about
 * @param lower the lower coordinate of the patch's side
 * @param upper the upper coordinate of the patch's side
 *
 * @return wether the boundary is a neumann boundary
 */
template <size_t D>
using IsNeumannFunc = std::function<bool(Side<D> s, const std::array<double, D> &lower,
                                         const std::array<double, D> &upper)>;

#endif
