#include "Cycle.h"
using namespace GMG;
void Cycle::prepCoarser(const Level &level)
{
	// calculate residual
	PW<Vec> r = level.getDomainCollection().getNewDomainVec();
	level.getOperator().apply(u_vectors.front(), r);
	VecAYPX(r, -1, f_vectors.front());
	// create vectors for coarser levels
	PW<Vec> new_u = level.getCoarser().getDomainCollection().getNewDomainVec();
	PW<Vec> new_f = level.getCoarser().getDomainCollection().getNewDomainVec();
	level.getRestrictor().restrict(new_f, r);
	u_vectors.push_front(new_u);
	f_vectors.push_front(new_f);
}
void Cycle::prepFiner(const Level &level)
{
	PW<Vec> old_u = u_vectors.front();
	u_vectors.pop_front();
	f_vectors.pop_front();
	level.getInterpolator().interpolate(old_u, u_vectors.front());
}
void Cycle::smooth(const Level &level)
{
	level.getSmoother().smooth(f_vectors.front(), u_vectors.front());
}
