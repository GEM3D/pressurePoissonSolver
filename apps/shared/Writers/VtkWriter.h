#ifndef VTKWRITER_H
#define VTKWRITER_H
#include "DomainCollection.h"
#include <map>
#include <set>
#include <string>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkMPIController.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLPMultiBlockDataWriter.h>
class VtkWriter
{
	private:
	DomainCollection<3>                                  dc;
	std::string                                       file_name;
	static vtkSmartPointer<vtkMultiProcessController> controller;
	std::map<int, vtkSmartPointer<vtkImageData>> images;
	std::set<vtkSmartPointer<vtkDoubleArray>>    arrays;
	vtkSmartPointer<vtkXMLPMultiBlockDataWriter> writer;
	vtkSmartPointer<vtkMultiBlockDataSet>        block;
	vtkSmartPointer<vtkMultiPieceDataSet>        data;

	public:
	VtkWriter(DomainCollection<3> &dc, std::string file_name);
	void add(Vec u, std::string name);
	void write();
};
#endif
