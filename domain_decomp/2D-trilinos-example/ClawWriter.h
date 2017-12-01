#ifndef CLAWWRITER_H
#define CLAWWRITER_H
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
class ClawWriter
{
	private:
	DomainSignatureCollection dsc;
	void writePatch(DomainSignature &d, std::ostream &os, vector_type &u, vector_type &resid);

	public:
	ClawWriter(DomainSignatureCollection &dsc);
	void write(vector_type &u, vector_type &resid);
};
#endif
