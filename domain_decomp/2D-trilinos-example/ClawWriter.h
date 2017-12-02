#ifndef CLAWWRITER_H
#define CLAWWRITER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
class ClawWriter
{
	private:
	DomainCollection dc;
	void writePatch(Domain &d, std::ostream &os, vector_type &u, vector_type &resid);

	public:
	ClawWriter(DomainCollection &dc);
	void write(vector_type &u, vector_type &resid);
};
#endif
