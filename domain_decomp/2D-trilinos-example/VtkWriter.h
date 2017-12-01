#ifndef VTKWRITER_H
#define VTKWRITER_H
#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
class VtkWriter
{
	private:
	DomainSignatureCollection dsc;

	public:
	VtkWriter(DomainSignatureCollection &dsc);
	void write(vector_type &u, vector_type &error, vector_type &resid);
};
#endif
