#ifndef CLAWWRITER_H
#define CLAWWRITER_H
#include "DomainCollection.h"
class ClawWriter
{
	private:
	DomainCollection<3> dc;
	void writePatch(Domain<3> &d, std::ostream &os, double *u_view, double *resid_view);

	public:
	ClawWriter(DomainCollection<3> &dc);
	void write(Vec u, Vec resid);
};
#endif
