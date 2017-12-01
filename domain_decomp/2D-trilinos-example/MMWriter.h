#ifndef MMWRITER_H
#define MMWRITER_H
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
#include <string>
class MMWriter
{
	private:
	DomainSignatureCollection dsc;
	bool                      amr;

	public:
	MMWriter(DomainSignatureCollection &dsc, bool amr);
	void write(const vector_type &u, std::string filename);
};
#endif
