#ifndef VTKWRITER_H
#define VTKWRITER_H
#include "DomainCollection.h"
#include "MyTypeDefs.h"
class VtkWriter
{
	private:
	DomainCollection dc;

	public:
	VtkWriter(DomainCollection &dc);
	void write(vector_type &u, vector_type &error, vector_type &resid);
};
#endif
