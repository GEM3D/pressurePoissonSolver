#ifndef MMWRITER_H
#define MMWRITER_H
#include "DomainCollection.h"
#include <string>
class MMWriter
{
	private:
	DomainCollection<3> dc;
	bool                amr;

	public:
	MMWriter(DomainCollection<3> &dc, bool amr);
	void write(const Vec u, std::string filename);
};
#endif
