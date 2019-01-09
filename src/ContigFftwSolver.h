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

#ifndef CONTIGFFTWSOLVER_H
#define CONTIGFFTWSOLVER_H
#include "DomainCollection.h"
#include <deque>
#include <fftw3.h>
struct FftwChunk {
	fftw_plan          plan;
	fftw_plan          plan_inv;
	std::deque<Domain> domains;
	int                n;
	int                size;
	void setup(double *eigs);
	void apply(double *fhat);
	void applyInv(double *fhat);
	void scale(double *fhat);
	virtual void get(double *fhat, double *tmp) = 0;
	virtual void put(double *fhat, double *tmp) = 0;
};
struct FftwChunkX : public FftwChunk {
	void get(double *fhat, double *tmp);
	void put(double *fhat, double *tmp);
};
struct FftwChunkY : public FftwChunk {
	void get(double *fhat, double *tmp);
	void put(double *fhat, double *tmp);
};
struct FftwChunkZ : public FftwChunk {
	void get(double *fhat, double *tmp);
	void put(double *fhat, double *tmp);
};
class ContigFftwSolver
{
	private:
        int num_cells;
	std::deque<FftwChunkX> x_chunks;
	std::deque<FftwChunkY> y_chunks;
	std::deque<FftwChunkZ> z_chunks;

	public:
    PW<Vec> eigs_vec;
	ContigFftwSolver() = default;
	ContigFftwSolver(DomainCollection &dc);
	void solve(Vec f , Vec u);
	void trans(Vec f , Vec u);
};
#endif
