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

#ifndef BICGSTABSOLVER_H
#define BICGSTABSOLVER_H
#include <Thunderegg/BiCGStab.h>
#include <Thunderegg/Domain.h>
#include <Thunderegg/PatchOperator.h>
#include <Thunderegg/PatchSolvers/PatchSolver.h>
#include <Thunderegg/ValVector.h>
#include <bitset>
#include <fftw3.h>
#include <map>
template <size_t D> class BiCGStabSolver : public PatchSolver<D>
{
	private:
	class SingleVG : public VectorGenerator<D>
	{
		private:
		std::array<int, D> lengths;

		public:
		SingleVG(const SchurInfo<D> &sinfo)
		{
			lengths = sinfo.pinfo->ns;
		}
		std::shared_ptr<Vector<D>> getNewVector()
		{
			return std::shared_ptr<Vector<D>>(new ValVector<D>(lengths));
		}
	};
	class SinglePatchVec : public Vector<D>
	{
		private:
		LocalData<D> ld;

		public:
		SinglePatchVec(const LocalData<D> &ld)
		{
			this->num_local_patches = 1;
			this->ld                = ld;
		}
		LocalData<D> getLocalData(int local_patch_id)
		{
			return ld;
		}
		const LocalData<D> getLocalData(int local_patch_id) const
		{
			return ld;
		}
	};
	class SinglePatchOp : public Operator<D>
	{
		private:
		std::shared_ptr<PatchOperator<D>> op;
		SchurInfo<D>                      sinfo;

		public:
		SinglePatchOp(const SchurInfo<D> &sinfo, std::shared_ptr<PatchOperator<D>> op)
		{
			this->sinfo = sinfo;
			this->op    = op;
		}
		void apply(std::shared_ptr<const Vector<D>> x, std::shared_ptr<Vector<D>> b) const
		{
			op->apply(sinfo, x->getLocalData(0), b->getLocalData(0));
		}
	};
	int                               n;
	bool                              initialized = false;
	static bool                       compareDomains();
	double                            lambda;
	std::array<int, D + 1>            npow;
	std::shared_ptr<PatchOperator<D>> op;
	int                               max_it;
	double                            tol;

	public:
	BiCGStabSolver(std::shared_ptr<PatchOperator<D>> op, double tol = 1e-12, int max_it = 1000)
	{
		this->op     = op;
		this->tol    = tol;
		this->max_it = max_it;
	}
	void solve(SchurInfo<D> &sinfo, std::shared_ptr<const Vector<D>> f,
	           std::shared_ptr<Vector<D>> u, std::shared_ptr<const Vector<D - 1>> gamma);
	void domainSolve(std::vector<SchurInfo<D>> &domains, std::shared_ptr<const Vector<D>> f,
	                 std::shared_ptr<Vector<D>> u, std::shared_ptr<const Vector<D - 1>> gamma)
	{
		for (SchurInfo<D> &sinfo : domains) {
			solve(sinfo, f, u, gamma);
		}
	}
	void addDomain(SchurInfo<D> &sinfo) {}
};
template <size_t D>
void BiCGStabSolver<D>::solve(SchurInfo<D> &sinfo, std::shared_ptr<const Vector<D>> f,
                              std::shared_ptr<Vector<D>>           u,
                              std::shared_ptr<const Vector<D - 1>> gamma)
{
	std::shared_ptr<SinglePatchOp>      single_op(new SinglePatchOp(sinfo, op));
	std::shared_ptr<VectorGenerator<D>> vg(new SingleVG(sinfo));

	std::shared_ptr<Vector<D>> f_single(
	new SinglePatchVec(f->getLocalData(sinfo.pinfo->local_index)));
	std::shared_ptr<Vector<D>> u_single(
	new SinglePatchVec(u->getLocalData(sinfo.pinfo->local_index)));

	auto f_copy = vg->getNewVector();
	f_copy->copy(f_single);
	op->addInterfaceToRHS(sinfo, gamma, f_copy->getLocalData(0));

	BiCGStab<D>::solve(vg, single_op, u_single, f_copy, nullptr, max_it, tol);
}
#endif
