#include "ClawWriter.h"
#include <fstream>
using namespace std;
ClawWriter::ClawWriter(DomainSignatureCollection &dsc) { this->dsc = dsc; }
void ClawWriter::write(vector_type &u, vector_type &resid)
{
	ofstream     t_file("fort.t0000");
	const string tab = "\t";
	t_file << 0.0 << tab << "time" << endl;
	t_file << 2 << tab << "meqn" << endl;
	t_file << dsc.domains.size() << tab << "ngrids" << endl;
	t_file << 2 << tab << "num_aux" << endl;
	t_file << 2 << tab << "num_dim" << endl;
	t_file.close();
	ofstream q_file("fort.q0000");

	q_file.precision(10);
	q_file << scientific;
	for (auto &p : dsc.domains) {
		DomainSignature &d = p.second;
		writePatch(d, q_file, u, resid);
	}
	q_file.close();
}
void ClawWriter::writePatch(DomainSignature &d, std::ostream &os, vector_type &u,
                            vector_type &resid)
{
	const string tab = "\t";
	os << d.id << tab << "grid_number" << endl;
	os << d.refine_level << tab << "AMR_level" << endl;
	os << 0 << tab << "block_number" << endl;
	os << 0 << tab << "mpi_rank" << endl;
	os << d.n << tab << "mx" << endl;
	os << d.n << tab << "my" << endl;
	os << d.x_start << tab << "xlow" << endl;
	os << d.y_start << tab << "ylow" << endl;
	int    n   = d.n;
	double h_x = d.x_length / n;
	double h_y = d.y_length / n;
	os << h_x << tab << "dx" << endl;
	os << h_y << tab << "dy" << endl;
	os << endl;
	int  start      = d.id * n * n;
	auto u_view     = u.get1dView();
	auto resid_view = resid.get1dView();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int loc = j + i * n;
			os << u_view[start + loc] << tab << resid_view[start + loc] * h_x * h_y << endl;
		}
		os << endl;
	}
}
