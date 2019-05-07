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

#ifndef SCHURMATRIXHELPER2D_H
#define SCHURMATRIXHELPER2D_H
#include <Thunderegg/SchurHelper.h>
#include <functional>
#include <memory>
#include <petscmat.h>
#include <valarray>
/**
 * @brief This class represents a collection of domains that a single processor owns.
 *
 * The purposes of this class:
 *   - Provide a member function for solving with a given interface vector.
 *   - Handle the initialization of the domains.
 *   - Provide member functions for calculating error, residual, etc.
 *   - Provide member functions that generate the Schur complement matrix.
 */
class SchurMatrixHelper2d
{
	private:
	std::shared_ptr<SchurHelper<2>> sh;
	int                             n;

	typedef std::function<void(int, int, std::shared_ptr<std::valarray<double>>, bool, bool)>
	     inserter;
	void assembleMatrix(inserter insertBlock);

	public:
	SchurMatrixHelper2d(std::shared_ptr<SchurHelper<2>> sh)
	{
		this->sh = sh;
		n        = sh->getLengths()[0];
	}
	/**
	 * @brief Form the Schur complement matrix
	 *
	 * @return the formed matrix
	 */
	PW_explicit<Mat> formCRSMatrix();
};
#endif
