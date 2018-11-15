#include "VtkWriter2d.h"
using namespace std;
vtkSmartPointer<vtkMultiProcessController> VtkWriter2d::controller;
VtkWriter2d::VtkWriter2d(DomainCollection<2> &dc, string file_name)
{
	if (controller == nullptr) {
		controller = vtkSmartPointer<vtkMPIController>::New();
		controller->Initialize(0, 0, 1);
	}
	this->dc        = dc;
	this->file_name = file_name;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	writer = vtkSmartPointer<vtkXMLPMultiBlockDataWriter>::New();
	writer->SetController(controller);
	block = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	data  = vtkSmartPointer<vtkMultiPieceDataSet>::New();
	data->SetNumberOfPieces(dc.num_global_domains);
	block->SetNumberOfBlocks(1);
	block->SetBlock(0, data);

	string out = file_name + ".vtmb";
	writer->SetFileName(out.c_str());
	writer->SetNumberOfPieces(1);
	writer->SetStartPiece(0);
	writer->SetInputData(block);
	if (rank == 0) { writer->SetWriteMetaFile(1); }
}
void VtkWriter2d::add(Vec u, string name)
{
	// create MultiPieceDataSet and fill with patch information

	double *u_view;
	VecGetArray(u, &u_view);

	for (auto &p : dc.domains) {
		Domain<2> &d     = *p.second;
		int        n     = d.n;
		double     h_x   = d.lengths[0] / n;
		double     h_y   = d.lengths[1] / n;
		int        start = d.id_local * n * n;

		// create image object
		vtkSmartPointer<vtkImageData> image = images[d.id];
		if (image == nullptr) {
			image        = vtkSmartPointer<vtkImageData>::New();
			images[d.id] = image;
			image->SetOrigin(d.starts[0], d.starts[1], 0);
			image->SetSpacing(h_x, h_y, 0);
			image->SetExtent(d.starts[0], d.starts[0] + d.lengths[0], d.starts[1],
			                 d.starts[1] + d.lengths[1], 0, 0);
			image->PrepareForNewData();
			image->SetDimensions(n + 1, n + 1, 1);
			// add image to dataset
			data->SetPiece(d.id_global, image);
		}

		// create solution vector
		vtkSmartPointer<vtkDoubleArray> solution = vtkSmartPointer<vtkDoubleArray>::New();
		solution->SetNumberOfComponents(1);
		solution->SetNumberOfValues(n * n);
		solution->SetName(name.c_str());
		for (int i = 0; i < n * n; i++) {
			solution->SetValue(i, u_view[start + i]);
		}
		image->GetCellData()->AddArray(solution);
	}

	VecRestoreArray(u, &u_view);
}
void VtkWriter2d::write()
{
	// write data
	writer->Update();
	writer->Write();
}
