#include "MMWriter.h"
#include <fstream>
using namespace std;
MMWriter::MMWriter(DomainCollection &dc, bool amr)
{
	this->amr = amr;
	this->dc  = dc;
}
void MMWriter::write(const Vec u, string filename)
{
	int      n = dc.n;
	ofstream os(filename);
	int      num_i, num_j, d_x;
	if (amr) {
		num_i = n * sqrt(dc.domains.size() / 5);
		num_j = n * sqrt(dc.domains.size() / 5);
		d_x   = sqrt(dc.domains.size() / 5);
	} else {
		num_i = n * sqrt(dc.domains.size());
		num_j = n * sqrt(dc.domains.size());
		d_x   = sqrt(dc.domains.size());
	}
	os << "%%MatrixMarket matrix array real general\n";
	os << num_i << ' ' << num_j << '\n';
	os.precision(15);
	double *u_view;
	VecGetArray(u, &u_view);
	for (int j = 0; j < num_j; j++) {
		int domain_j   = j / n;
		int internal_j = j % n;
		for (int i = 0; i < num_i; i++) {
			int domain_i   = i / n;
			int internal_i = i % n;
			int id         = domain_i * d_x + domain_j;
			int start      = n * n * dc.domains[id].id_local;
			os << u_view[start + internal_i * n + internal_j] << '\n';
		}
	}
	os.close();
	if (amr) {
		ofstream os(filename + ".amr");
		int      num_i = 2 * n * sqrt(dc.domains.size() / 5);
		int      num_j = 2 * n * sqrt(dc.domains.size() / 5);
		int      d_x   = 2 * sqrt(dc.domains.size() / 5);
		os << "%%MatrixMarket matrix array real general\n";
		os << num_i << ' ' << num_j << '\n';
		os.precision(15);
		for (int j = 0; j < num_j; j++) {
			int domain_j   = j / n;
			int internal_j = j % n;
			for (int i = 0; i < num_i; i++) {
				int domain_i   = i / n;
				int internal_i = i % n;
				int id         = d_x * d_x / 4 + domain_i * d_x + domain_j;
				int start      = n * n * dc.domains[id].id_local;
				os << u_view[start + internal_i * n + internal_j] << '\n';
			}
		}
		os.close();
	}
	VecRestoreArray(u, &u_view);
}
