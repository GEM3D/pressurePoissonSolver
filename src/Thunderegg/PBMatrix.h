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

#ifndef PBMATRIX_H
#define PBMATRIX_H
#include "PW.h"
#include <functional>
#include <map>
#include <memory>
#include <petscmat.h>
#include <petscpc.h>
#include <set>
#include <valarray>
struct PBlock {
	int                                    i;
	int                                    j;
	std::shared_ptr<std::valarray<double>> coeffs;
	std::function<int(int, int, int)> col_trans;
	std::function<int(int, int, int)> row_trans;
	friend bool operator<(const PBlock &l, const PBlock &r)
	{
		char lc = l.col_trans(2, 0, 0);
		char lr = l.row_trans(2, 0, 1);
		char rc = r.col_trans(2, 0, 0);
		char rr = r.row_trans(2, 0, 1);
		return std::tie(l.coeffs, l.i, l.j, lc, lr) < std::tie(r.coeffs, r.i, r.j, rc, rr);
	}
};
struct BlockJacobiSmoother;
class PBMatrix
{
	private:
	int              n;
	std::set<PBlock> blocks;
	typedef std::function<int(int, int, int)> trans_type;
	typedef std::map<trans_type, std::valarray<int>, std::function<bool(trans_type, trans_type)>>
	         map_type;
	map_type trans;
	int      local_size  = 0;
	int      global_size = 0;

	public:
	PBMatrix(int n, int local_size, int global_size);
	void insertBlock(int i, int j, std::shared_ptr<std::valarray<double>> coeffs,
	                 std::function<int(int, int, int)>                    col_trans,
	                 std::function<int(int, int, int)>                    row_trans);
	void apply(Vec x, Vec b);
	void       finalize();
	PBMatrix * getDiagInv();
	BlockJacobiSmoother getBlockJacobiSmoother();
	static int multiply(Mat A, Vec x, Vec b)
	{
		PBMatrix *pba = nullptr;
		MatShellGetContext(A, &pba);
		pba->apply(x, b);
		return 0;
	}
	static int destroy(Mat A)
	{
		PBMatrix *pba = nullptr;
		MatShellGetContext(A, &pba);
		delete pba;
		return 0;
	}
	PW_explicit<Mat> getMatrix()
	{
		PW<Mat> A;
		MatCreateShell(MPI_COMM_WORLD, local_size, local_size, global_size, global_size, this, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) multiply);
		MatShellSetOperation(A, MATOP_DESTROY, (void (*)(void)) destroy);
		return A;
	}
	static int multiplyPC(PC A, Vec x, Vec b)
	{
		PBMatrix *pba = nullptr;
		PCShellGetContext(A, (void **) &pba);
		pba->apply(x, b);
		return 0;
	}
	void getPrec(PC P)
	{
		PCSetType(P, PCSHELL);
		PCShellSetContext(P, this);
		PCShellSetApply(P, multiplyPC);
	}
};
struct BlockJacobiSmoother {
	std::shared_ptr<PBMatrix> D;
	std::shared_ptr<PBMatrix> R;
	void apply(Vec x, Vec b)
	{
		PW<Vec> tmp;
		VecDuplicate(x, &tmp);
		R->apply(x, tmp);
		VecAYPX(tmp, -1, b);
		D->apply(tmp, x);
	}
};
#endif
