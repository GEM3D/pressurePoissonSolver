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
