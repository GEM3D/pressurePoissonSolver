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

#include "ClawWriter.h"
#include <fstream>
using namespace std;
ClawWriter::ClawWriter(std::shared_ptr<Domain<2>> domain)
{
	this->domain = domain;
}
void ClawWriter::write(Vec u, Vec resid)
{
	ofstream     t_file("fort.t0000");
	const string tab = "\t";
	t_file << 0.0 << tab << "time" << endl;
	t_file << 2 << tab << "meqn" << endl;
	t_file << domain->getNumLocalPatches() << tab << "ngrids" << endl;
	t_file << 2 << tab << "num_aux" << endl;
	t_file << 2 << tab << "num_dim" << endl;
	t_file.close();
	ofstream q_file("fort.q0000");

	double *u_view, *resid_view;
	VecGetArray(u, &u_view);
	VecGetArray(resid, &resid_view);
	q_file.precision(10);
	q_file << scientific;
	for (auto &pinfo : domain->getPatchInfoVector()) {
		writePatch(*pinfo, q_file, u_view, resid_view);
	}
	VecRestoreArray(u, &u_view);
	VecRestoreArray(resid, &resid_view);
	q_file.close();
}
void ClawWriter::writePatch(PatchInfo<2> &d, std::ostream &os, double *u_view, double *resid_view)
{
	const string tab = "\t";
	os << d.id << tab << "grid_number" << endl;
	os << d.refine_level << tab << "AMR_level" << endl;
	os << 0 << tab << "block_number" << endl;
	os << 0 << tab << "mpi_rank" << endl;
	os << d.ns[0] << tab << "mx" << endl;
	os << d.ns[1] << tab << "my" << endl;
	os << d.starts[0] << tab << "xlow" << endl;
	os << d.starts[1] << tab << "ylow" << endl;
	os << d.spacings[0] << tab << "dx" << endl;
	os << d.spacings[1] << tab << "dy" << endl;
	os << endl;
	int start = d.id * d.ns[0] * d.ns[1];
	for (int i = 0; i < d.ns[0]; i++) {
		for (int j = 0; j < d.ns[1]; j++) {
			int loc = j + i * d.ns[1];
			os << u_view[start + loc] << tab
			   << resid_view[start + loc] * d.spacings[0] * d.spacings[1] << endl;
		}
		os << endl;
	}
}
