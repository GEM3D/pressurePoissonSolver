#ifndef CLAWWRITER_H
#define CLAWWRITER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
class ClawWriter
{
	private:
	DomainCollection dc;
	void writePatch(Domain &d, std::ostream &os, double *u_view, double *resid_view);

	public:
	ClawWriter(DomainCollection &dc);
	void write(Vec u, Vec resid);
};
#endif
