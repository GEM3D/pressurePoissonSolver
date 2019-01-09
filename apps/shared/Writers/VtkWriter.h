/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively 
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

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
