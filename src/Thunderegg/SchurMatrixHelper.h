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

#include <Thunderegg/Experimental/PBMatrix.h>
#include <Thunderegg/SchurHelper.h>
#include <functional>
#include <petscmat.h>
#include <valarray>
struct Block;
class SchurMatrixHelper
{
	private:
	std::shared_ptr<SchurHelper<3>> sh;
	int                             n;

	typedef std::function<void(Block *, std::shared_ptr<std::valarray<double>>)> inserter;
	void assembleMatrix(inserter insertBlock);

	public:
	SchurMatrixHelper(std::shared_ptr<SchurHelper<3>> sh)
	{
		this->sh = sh;
        n=sh->getLengths()[0];
	}
	PW_explicit<Mat> formCRSMatrix();
	PBMatrix *       formPBMatrix();
	void             getPBDiagInv(PC p);
	PW_explicit<Mat> getPBMatrix();
	PW_explicit<Mat> getPBDiagInv();
};
