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

#ifndef PETSCVECTOR_H
#define PETSCVECTOR_H
#include <Vector.h>
#include <petscvec.h>
class PetscLDM : public LocalDataManager
{
	private:
	// PW<Vec> vec;
	Vec     vec;
	double *vec_view;
	bool    is_read_only;

	public:
	PetscLDM(Vec vec, bool is_read_only)
	{
		this->vec          = vec;
		this->is_read_only = is_read_only;
		if (is_read_only) {
			VecGetArrayRead(vec, const_cast<const double **>(&vec_view));
		} else {
			VecGetArray(vec, &vec_view);
		}
	}
	~PetscLDM()
	{
		if (is_read_only) {
			VecRestoreArrayRead(vec, const_cast<const double **>(&vec_view));
		} else {
			VecRestoreArray(vec, &vec_view);
		}
	}
	double *getVecView() const
	{
		return vec_view;
	}
};
template <size_t D> class PetscVector : public Vector<D>

{
	private:
	int                patch_stride;
	bool               own;
	std::array<int, D> strides;
	std::array<int, D> lengths;

	public:
	Vec vec;
	PetscVector(Vec vec, const std::array<int, D> &lengths, bool own = true)
	{
		this->own     = own;
		this->lengths = lengths;
		this->vec     = vec;
		strides[0]    = 1;
		for (size_t i = 1; i < D; i++) {
			strides[i] = lengths[i - 1] * strides[i - 1];
		}
		patch_stride = strides[D - 1] * lengths[D - 1];
		VecGetLocalSize(vec, &this->num_local_patches);
		this->num_local_patches /= patch_stride;
	}
	~PetscVector()
	{
		if (own) { VecDestroy(&vec); }
	}
	LocalData<D> getLocalData(int local_patch_id)
	{
		std::shared_ptr<PetscLDM> ldm(new PetscLDM(vec, true));
		double *                  data = ldm->getVecView() + patch_stride * local_patch_id;
		return LocalData<D>(data, strides, lengths, ldm);
	}
	const LocalData<D> getLocalData(int local_patch_id) const
	{
		std::shared_ptr<PetscLDM> ldm(new PetscLDM(vec, false));
		double *                  data = ldm->getVecView() + patch_stride * local_patch_id;
		return LocalData<D>(data, strides, lengths, std::move(ldm));
	}
	void set(double alpha)
	{
		VecSet(vec, alpha);
	}
	void scale(double alpha)
	{
		VecScale(vec, alpha);
	}
	void shift(double delta)
	{
		VecShift(vec, delta);
	}
    void add(std::shared_ptr<const Vector<D>> b)
	{
		const PetscVector<D> *b_vec = dynamic_cast<const PetscVector<D> *>(b.get());
		if (b_vec != nullptr) {
			VecAXPY(vec, 1, b_vec->vec);
		} else {
			Vector<D>::add(b);
		}
	}
	void scaleThenAdd(double alpha, std::shared_ptr<const Vector<D>> b)
	{
		const PetscVector<D> *b_vec = dynamic_cast<const PetscVector<D> *>(b.get());
		if (b_vec != nullptr) {
			VecAYPX(vec, alpha, b_vec->vec);
		} else {
			Vector<D>::scaleThenAdd(alpha, b);
		}
	}
};
#endif
