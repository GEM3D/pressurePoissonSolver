#include "GMGHelper.h"
#include "GMGAvgRstr.h"
#include "GMGDrctIntp.h"
using namespace std;
GMGHelper::GMGHelper(int n, OctTree t, DomainCollection &dc, SchurHelper &sh)
{
	num_levels = t.num_levels;
	top_level  = t.num_levels - 1;
	// generate and balance levels
	levels.resize(num_levels);
	levels[top_level] = dc;
	for (int i = 0; i < top_level; i++) {
		levels[i]   = DomainCollection(t, i + 1);
		levels[i].n = n;
		for (auto &p : levels[i].domains) {
			p.second.n = n;
		}
	}
	// rebalance
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size > 1) {
		for (int i = top_level - 1; i >= 0; i--) {
			levels[i].zoltanBalanceWithLower(levels[i + 1]);
		}
	}
	shs.resize(num_levels);
	u_vectors.resize(num_levels);
	f_vectors.resize(num_levels);
	r_vectors.resize(num_levels);
	restrictors.resize(num_levels - 1);
	interpolators.resize(num_levels - 1);
	shs[top_level]       = sh;
	r_vectors[top_level] = dc.getNewDomainVec();
	for (int i = 0; i < top_level; i++) {
		restrictors[i]   = new GMGAvgRstr(levels[i], levels[i + 1]);
		interpolators[i] = new GMGDrctIntp(levels[i], levels[i + 1]);
		shs[i]           = SchurHelper(levels[i], sh.getSolver(), sh.getOp(), sh.getInterpolator());
		f_vectors[i]     = levels[i].getNewDomainVec();
		u_vectors[i]     = levels[i].getNewDomainVec();
		r_vectors[i]     = levels[i].getNewDomainVec();
	}
}
void GMGHelper::apply(Vec f, Vec u)
{
	f_vectors[top_level] = f;
	u_vectors[top_level] = u;
	// finest level
	// smooth
	shs[top_level].solveWithSolution(f, u);
	// residual
	shs[top_level].apply(u, r_vectors[top_level]);
	VecAYPX(r_vectors[top_level], -1, f);
	// down cycle
	for (int i = top_level - 1; i >= 1; i--) {
		restrictors[i]->restrict(f_vectors[i], r_vectors[i + 1]);
		// smooth
		VecScale(u_vectors[i], 0);
		shs[i].solveWithSolution(f_vectors[i], u_vectors[i]);
		// residual
		shs[i].apply(u_vectors[i], r_vectors[i]);
		VecAYPX(r_vectors[i], -1, f_vectors[i]);
	}
	// coarse level
	restrictors[0]->restrict(f_vectors[0], r_vectors[1]);
	// smooth
	VecScale(u_vectors[0], 0);
	shs[0].solveWithSolution(f_vectors[0], u_vectors[0]);
	interpolators[0]->interpolate(u_vectors[0], u_vectors[1]);
	// up cycle
	for (int i = 1; i <= num_levels - 2; i++) {
		// smooth
		shs[i].solveWithSolution(f_vectors[i], u_vectors[i]);
	    interpolators[i]->interpolate(u_vectors[i], u_vectors[i+1]);
	}
	// finest level
	// smooth
	shs[top_level].solveWithSolution(f, u);
}
