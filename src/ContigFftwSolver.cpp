#include "ContigFftwSolver.h"
#include <set>
#include <valarray>
using namespace std;
void FftwChunk::setup(double *eigs)
{
	n    = domains.front().n;
	size = domains.size();
	valarray<double> in(size * n * n * n);
	valarray<double> out(size * n * n * n);
	fftw_r2r_kind    transform     = FFTW_RODFT10;
	fftw_r2r_kind    transform_inv = FFTW_RODFT01;
	int              num           = n * size;

	plan     = fftw_plan_many_r2r(1, &num, n * n, &in[0], nullptr, 1, n * size, &out[0], nullptr, 1,
                              n * size, &transform, FFTW_MEASURE | FFTW_DESTROY_INPUT);
	plan_inv = fftw_plan_many_r2r(1, &num, n * n, &in[0], nullptr, 1, n * size, &out[0], nullptr, 1,
	                              n * size, &transform_inv, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	double h = domains.front().x_length / n;
	get(eigs, &in[0]);
	valarray<double> ones(n * n);
	ones = 1;
	for (int xi = 0; xi < size * n; xi++) {
		in[slice(xi, n * n, size * n)]
		-= 4 / (h * h) * pow(sin((xi + 1) * M_PI / (2 * size * n)), 2) * ones;
	}
	put(eigs, &in[0]);
}
void FftwChunk::apply(double *fhat)
{
	valarray<double> in(size * n * n * n);
	valarray<double> out(size * n * n * n);
	get(fhat, &in[0]);
	// in /= size;
	in /= sqrt(2 * size * n);
	fftw_execute_r2r(plan, &in[0], &out[0]);
	put(fhat, &out[0]);
}
void FftwChunk::applyInv(double *fhat)
{
	valarray<double> in(size * n * n * n);
	valarray<double> out(size * n * n * n);
	get(fhat, &in[0]);
	// in /= size;
	// in /= sqrt(2 * size * n);
	in /= sqrt(2 * size * n);
	fftw_execute_r2r(plan_inv, &in[0], &out[0]);
	// out /= sqrt(2 * size * n);
	put(fhat, &out[0]);
}
void FftwChunk::scale(double *fhat)
{
	valarray<double> in(size * n * n * n);
	get(fhat, &in[0]);
	// in /= 8 * size * n*n*n;
	put(fhat, &in[0]);
}
void FftwChunkX::get(double *fhat, double *tmp)
{
	int pos = 0;
	for (Domain d : domains) {
		int start_idx = d.id_local * n * n * n;
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					tmp[xi + pos * n + yi * n * size + zi * n * n * size]
					= fhat[start_idx + xi + yi * n + zi * n * n];
				}
			}
		}
		pos++;
	}
}
void FftwChunkX::put(double *fhat, double *tmp)
{
	int pos = 0;
	for (Domain d : domains) {
		int start_idx = d.id_local * n * n * n;
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					fhat[start_idx + xi + yi * n + zi * n * n]
					= tmp[xi + pos * n + yi * n * size + zi * n * n * size];
				}
			}
		}
		pos++;
	}
}
void FftwChunkY::get(double *fhat, double *tmp)
{
	int pos = 0;
	for (Domain d : domains) {
		int start_idx = d.id_local * n * n * n;
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					tmp[yi + pos * n + zi * n * size + xi * n * n * size]
					= fhat[start_idx + xi + yi * n + zi * n * n];
				}
			}
		}
		pos++;
	}
}
void FftwChunkY::put(double *fhat, double *tmp)
{
	int pos = 0;
	for (Domain d : domains) {
		int start_idx = d.id_local * n * n * n;
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					fhat[start_idx + xi + yi * n + zi * n * n]
					= tmp[yi + pos * n + zi * n * size + xi * n * n * size];
				}
			}
		}
		pos++;
	}
}
void FftwChunkZ::get(double *fhat, double *tmp)
{
	int pos = 0;
	for (Domain d : domains) {
		int start_idx = d.id_local * n * n * n;
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					tmp[zi + pos * n + xi * n * size + yi * n * n * size]
					= fhat[start_idx + xi + yi * n + zi * n * n];
				}
			}
		}
		pos++;
	}
}
void FftwChunkZ::put(double *fhat, double *tmp)
{
	int pos = 0;
	for (Domain d : domains) {
		int start_idx = d.id_local * n * n * n;
		for (int zi = 0; zi < n; zi++) {
			for (int yi = 0; yi < n; yi++) {
				for (int xi = 0; xi < n; xi++) {
					fhat[start_idx + xi + yi * n + zi * n * n]
					= tmp[zi + pos * n + xi * n * size + yi * n * n * size];
				}
			}
		}
		pos++;
	}
}
ContigFftwSolver::ContigFftwSolver(DomainCollection &dc)
{
	eigs_vec = dc.getNewDomainVec();
	double *eigs;
	VecGetArray(eigs_vec, &eigs);
	set<shared_ptr<Domain>> todo;
	for (auto p : dc.domains) {
		todo.insert(p.second);
	}
	auto getNbrId = [&](int id, Side s) {
		Domain &d      = *dc.domains[id];
		int     retval = -1;
		if (d.hasNbr(s)) { retval = d.getNormalNbrInfo(s).id; }
		return retval;
	};
	while (!todo.empty()) {
		Domain     d = **todo.begin();
		FftwChunkX chunk;
		// east nbrs
		for (int id = d.id; id != -1; id = getNbrId(id, Side::east)) {
			shared_ptr<Domain> nbr = dc.domains[id];
			chunk.domains.push_back(*nbr);
			todo.erase(nbr);
		}
		// west nbrs
		if (d.hasNbr(Side::west)) {
			for (int id = getNbrId(d.id, Side::west); id != -1; id = getNbrId(id, Side::west)) {
				shared_ptr<Domain> nbr = dc.domains[id];
				chunk.domains.push_front(*nbr);
				todo.erase(nbr);
			}
		}
		chunk.setup(eigs);
		x_chunks.push_back(chunk);
	}
	for (auto p : dc.domains) {
		todo.insert(p.second);
	}
	while (!todo.empty()) {
		Domain     d = **todo.begin();
		FftwChunkY chunk;
		// north nbrs
		for (int id = d.id; id != -1; id = getNbrId(id, Side::north)) {
			shared_ptr<Domain> nbr = dc.domains[id];
			chunk.domains.push_back(*nbr);
			todo.erase(nbr);
		}
		// south nbrs
		if (d.hasNbr(Side::south)) {
			for (int id = getNbrId(d.id, Side::south); id != -1; id = getNbrId(id, Side::south)) {
				shared_ptr<Domain> nbr = dc.domains[id];
				chunk.domains.push_front(*nbr);
				todo.erase(nbr);
			}
		}
		chunk.setup(eigs);
		y_chunks.push_back(chunk);
	}
	for (auto p : dc.domains) {
		todo.insert(p.second);
	}
	while (!todo.empty()) {
		Domain     d = **todo.begin();
		FftwChunkZ chunk;
		chunk.n = d.n;
		// top nbrs
		for (int id = d.id; id != -1; id = getNbrId(id, Side::top)) {
			shared_ptr<Domain> nbr = dc.domains[id];
			chunk.domains.push_back(*nbr);
			todo.erase(nbr);
		}
		// bottom nbrs
		if (d.hasNbr(Side::bottom)) {
			for (int id = getNbrId(d.id, Side::bottom); id != -1; id = getNbrId(id, Side::bottom)) {
				shared_ptr<Domain> nbr = dc.domains[id];
				chunk.domains.push_front(*nbr);
				todo.erase(nbr);
			}
		}
		chunk.setup(eigs);
		z_chunks.push_back(chunk);
	}
	VecRestoreArray(eigs_vec, &eigs);
	int n     = dc.getN();
	num_cells = dc.domains.size() * n * n * n;
}
void ContigFftwSolver::solve(Vec f, Vec u)
{
	VecAYPX(u, 0, f);
	double *fhat;
	VecGetArray(u, &fhat);
	for (FftwChunkX c : x_chunks) {
		c.scale(fhat);
		c.apply(fhat);
	}
	for (FftwChunkY c : y_chunks) {
		c.apply(fhat);
	}
	for (FftwChunkZ c : z_chunks) {
		c.apply(fhat);
	}
	VecPointwiseDivide(u, u, eigs_vec);
	for (FftwChunkZ c : z_chunks) {
		c.applyInv(fhat);
	}
	for (FftwChunkY c : y_chunks) {
		c.applyInv(fhat);
	}
	for (FftwChunkX c : x_chunks) {
		c.applyInv(fhat);
	}
	VecRestoreArray(u, &fhat);
}
void ContigFftwSolver::trans(Vec f, Vec u)
{
	VecAYPX(u, 0, f);
	double *fhat;
	VecGetArray(u, &fhat);
	for (FftwChunkX c : x_chunks) {
		c.apply(fhat);
	}
	for (FftwChunkY c : y_chunks) {
		c.apply(fhat);
	}
	for (FftwChunkZ c : z_chunks) {
		c.apply(fhat);
	}
	VecRestoreArray(u, &fhat);
}
