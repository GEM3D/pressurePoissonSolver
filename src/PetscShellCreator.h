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

#ifndef PETSCPCSHELLOP_H
#define PETSCPCSHELLOP_H
#include <Operators/Operator.h>
/**
 * @brief Base class for operators
 */
class PetscShellCreator
{
	private:
	template <size_t D> class PetscPCShellOpDomain
	{
		private:
		std::shared_ptr<Operator<D>>         op;
		std::shared_ptr<DomainCollection<D>> dc;

		public:
		PetscPCShellOpDomain(std::shared_ptr<Operator<D>>         op,
		                     std::shared_ptr<DomainCollection<D>> dc)
		{
			this->op = op;
			this->dc = dc;
		}
		static int applyMat(Mat A, Vec x, Vec b)
		{
			PetscPCShellOpDomain<D> *wrap = nullptr;
			MatShellGetContext(A, &wrap);
			std::shared_ptr<Vector<D>> x_vec(new PetscVector<D>(x, wrap->dc->getLengths(), false));
			std::shared_ptr<Vector<D>> b_vec(new PetscVector<D>(b, wrap->dc->getLengths(), false));
			wrap->op->apply(x_vec, b_vec);
			return 0;
		}
		static int destroyMat(Mat A)
		{
			PetscPCShellOpDomain<D> *wrap = nullptr;
			MatShellGetContext(A, &wrap);
			delete wrap;
			return 0;
		}
		static int apply(PC P, Vec x, Vec b)
		{
			PetscPCShellOpDomain<D> *wrap = nullptr;
			PCShellGetContext(P, (void **) &wrap);
			std::shared_ptr<Vector<D>> x_vec(new PetscVector<D>(x, wrap->dc->getLengths(), false));
			std::shared_ptr<Vector<D>> b_vec(new PetscVector<D>(b, wrap->dc->getLengths(), false));
			wrap->op->apply(x_vec, b_vec);
			return 0;
		}
		static int destroy(PC P)
		{
			PetscPCShellOpDomain<D> *wrap = nullptr;
			PCShellGetContext(P, (void **) &wrap);
			delete wrap;
			return 0;
		}
	};
	template <size_t D> class PetscPCShellOpSchur
	{
		private:
		std::shared_ptr<Operator<D>>        op;
		std::shared_ptr<SchurHelper<D + 1>> sh;

		public:
		PetscPCShellOpSchur(std::shared_ptr<Operator<D>> op, std::shared_ptr<SchurHelper<D + 1>> sh)
		{
			this->op = op;
			this->sh = sh;
		}
		static int applyMat(Mat A, Vec x, Vec b)
		{
			PetscPCShellOpSchur<D> *wrap = nullptr;
			MatShellGetContext(A, &wrap);
			std::shared_ptr<Vector<D>> x_vec(new PetscVector<D>(x, wrap->sh->getLengths(), false));
			std::shared_ptr<Vector<D>> b_vec(new PetscVector<D>(b, wrap->sh->getLengths(), false));
			wrap->op->apply(x_vec, b_vec);
			return 0;
		}
		static int destroyMat(Mat A)
		{
			PetscPCShellOpSchur<D> *wrap = nullptr;
			MatShellGetContext(A, &wrap);
			delete wrap;
			return 0;
		}
		static int apply(PC P, Vec x, Vec b)
		{
			PetscPCShellOpSchur<D> *wrap = nullptr;
			PCShellGetContext(P, (void **) &wrap);
			std::shared_ptr<Vector<D>> x_vec(new PetscVector<D>(x, wrap->sh->getLengths(), false));
			std::shared_ptr<Vector<D>> b_vec(new PetscVector<D>(b, wrap->sh->getLengths(), false));
			wrap->op->apply(x_vec, b_vec);
			return 0;
		}
		static int destroy(PC P)
		{
			PetscPCShellOpSchur<D> *wrap = nullptr;
			PCShellGetContext(P, (void **) &wrap);
			delete wrap;
			return 0;
		}
	};

	public:
	template <size_t D>
	static void getPCShell(PC pc, std::shared_ptr<Operator<D>> op,
	                       std::shared_ptr<DomainCollection<D>> dc)
	{
		PetscPCShellOpDomain<D> *wrap = new PetscPCShellOpDomain<D>(op, dc);
		PCSetType(pc, PCSHELL);
		PCShellSetContext(pc, wrap);
		PCShellSetApply(pc, PetscPCShellOpDomain<D>::apply);
		PCShellSetDestroy(pc, PetscPCShellOpDomain<D>::destroy);
	}
	template <size_t D>
	static void getPCShell(PC pc, std::shared_ptr<Operator<D>> op,
	                       std::shared_ptr<SchurHelper<D + 1>> sh)
	{
		PetscPCShellOpSchur<D> *wrap = new PetscPCShellOpSchur<D>(op, sh);
		PCSetType(pc, PCSHELL);
		PCShellSetContext(pc, wrap);
		PCShellSetApply(pc, PetscPCShellOpSchur<D>::apply);
		PCShellSetDestroy(pc, PetscPCShellOpSchur<D>::destroy);
	}
	template <size_t D>
	static PW_explicit<Mat> getMatShell(std::shared_ptr<Operator<D>>         op,
	                                    std::shared_ptr<DomainCollection<D>> dc)
	{
		PetscPCShellOpDomain<D> *wrap = new PetscPCShellOpDomain<D>(op, dc);
		int                      M    = dc->getGlobalNumCells();
		int                      m    = dc->getLocalNumCells();
		PW<Mat>                  A;
		MatCreateShell(MPI_COMM_WORLD, m, m, M, M, wrap, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) PetscPCShellOpDomain<D>::applyMat);
		MatShellSetOperation(A, MATOP_DESTROY,
		                     (void (*)(void)) PetscPCShellOpDomain<D>::destroyMat);
		return A;
	}
	template <size_t D>
	static PW_explicit<Mat> getMatShell(std::shared_ptr<Operator<D>>    op,
	                                    std::shared_ptr<SchurHelper<D>> sh)
	{
		PetscPCShellOpSchur<D> *wrap = new PetscPCShellOpSchur<D>(op, sh);
		int                     M    = sh->getSchurVecGlobalSize();
		int                     m    = sh->getSchurVecLocalSize();
		PW<Mat>                 A;
		MatCreateShell(MPI_COMM_WORLD, m, m, M, M, wrap, &A);
		MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) PetscPCShellOpSchur<D>::applyMat);
		MatShellSetOperation(A, MATOP_DESTROY, (void (*)(void)) PetscPCShellOpSchur<D>::destroyMat);
		return A;
	}
};
#endif
