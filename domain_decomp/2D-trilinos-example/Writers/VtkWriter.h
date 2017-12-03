#ifndef VTKWRITER_H
#define VTKWRITER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
#include <string>
class VtkWriter
{
	private:
	DomainCollection dc;

	public:
	VtkWriter(DomainCollection &dc);
	void write(std::string file_name, vector_type &u, vector_type &error, vector_type &resid);
};
#endif
