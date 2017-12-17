#ifndef MMWRITER_H
#define MMWRITER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
#include <string>
class MMWriter
{
	private:
	DomainCollection dc;
	bool             amr;

	public:
	MMWriter(DomainCollection &dc, bool amr);
	void write(const Vec u, std::string filename);
};
#endif
