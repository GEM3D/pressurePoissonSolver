#include "DomainSignatureCollection.h"
#include <iostream>
#include <fstream>
#include <mpi.h>
using namespace std;
DomainSignatureCollection::DomainSignatureCollection(string file_name, int rank){
	this->rank = rank;
	if (rank == 0) {
		num_global_interfaces = 0;
		ifstream mesh(file_name);
		while (!mesh.eof()) {
			DomainSignature ds;
			int             null, nl, nr, el, er, sl, sr, wl, wr;
			mesh >> ds.id >> null >> wl >> wr >> er >> el >> sr >> sl >> nl >> nr;

			// north
			ds.nbr(Side::north) = nl;
			if (nl != nr) {
				ds.nbrRight(Side::north) = nr;
				ds.setHasFineNbr(Side::north);
			}

			// east
			ds.nbr(Side::east) = el;
			if (el != er) {
				ds.nbrRight(Side::east) = er;
				ds.setHasFineNbr(Side::east);
			}

			// south
			ds.nbr(Side::south) = sl;
			if (sl != sr) {
				ds.nbrRight(Side::south) = sr;
				ds.setHasFineNbr(Side::south);
			}

			// west
			ds.nbr(Side::west) = wl;
			if (wl != wr) {
				ds.nbrRight(Side::west) = wr;
				ds.setHasFineNbr(Side::west);
			}
			domains[ds.id] = ds;
		}
		determineCoarseness();
		determineAmrLevel();
		determineXY();
		num_global_domains = domains.size();
	}
	MPI_Bcast(&num_global_domains, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
void DomainSignatureCollection::determineCoarseness(){
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	int        first = domains.begin()->first;
	queue.push_back(first);
	enqueued.insert(first);
	while (!queue.empty()) {
		int              curr = queue.front();
		DomainSignature &d    = domains.at(curr);
		queue.pop_front();
		visited.insert(curr);
        Side s = Side::north;
        do{
			if (d.hasNbr(s)) {
				if (d.hasFineNbr(s)) {
					DomainSignature &nbr_left  = domains.at(d.nbr(s));
					DomainSignature &nbr_right = domains.at(d.nbrRight(s));

					// set center indexes
					nbr_left.setHasCoarseNbr(!s);
					nbr_left.setLeftOfCoarse(!s);
					nbr_right.setHasCoarseNbr(!s);

					// enqueue domains
					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
					if (enqueued.count(d.nbrRight(s)) == 0) {
						queue.push_back(d.nbrRight(s));
						enqueued.insert(d.nbrRight(s));
					}
				} else {
					// enqueue domain
					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
				}
			}
			s++;
		} while (s != Side::north);
	}
}
void DomainSignatureCollection::determineAmrLevel(){
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	int        first = domains.begin()->first;
	queue.push_back(first);
	enqueued.insert(first);
	int min_level = 1;
	while (!queue.empty()) {
		int              curr = queue.front();
		DomainSignature &d    = domains.at(curr);
		int              curr_level = d.refine_level;
		queue.pop_front();
		visited.insert(curr);
        Side s = Side::north;
        do{
			if (d.hasNbr(s) && visited.count(d.nbr(s)) == 0) {
				// fine case
				if (d.hasFineNbr(s)) {
					DomainSignature &nbr_left  = domains.at(d.nbr(s));
					DomainSignature &nbr_right = domains.at(d.nbrRight(s));

					nbr_left.refine_level      = curr_level + 1;
					nbr_right.refine_level = curr_level + 1;
					// enqueue domains
					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
					if (enqueued.count(d.nbrRight(s)) == 0) {
						queue.push_back(d.nbrRight(s));
						enqueued.insert(d.nbrRight(s));
					}
					// coarse case
				} else if (d.hasCoarseNbr(s)) {
					DomainSignature &nbr = domains.at(d.nbr(s));

					nbr.refine_level     = curr_level - 1;
					if (curr_level - 1 < min_level) {
						min_level = curr_level;
					}

					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
					// normal case
				} else {
					DomainSignature &nbr = domains.at(d.nbr(s));
					nbr.refine_level     = curr_level;
					// enqueue domain
					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
				}
			}
			s++;
		} while (s != Side::north);
	}
}
void DomainSignatureCollection::determineXY(){
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	int        first = domains.begin()->first;
	queue.push_back(first);
	enqueued.insert(first);
	double x_min = 0;
	double y_min = 0;
	double x_max = 0;
	double y_max = 0;
	while (!queue.empty()) {
		int              curr = queue.front();
		DomainSignature &d    = domains.at(curr);
		queue.pop_front();
		if (d.x_start < x_min) {
			x_min = d.x_start;
		}
		if (d.y_start < y_min) {
			y_min = d.y_start;
		}
        if(d.x_start+d.x_length>x_max){
            x_max = d.x_start+d.x_length;
        }
        if(d.y_start+d.y_length>y_max){
            y_max = d.y_start+d.y_length;
        }
		visited.insert(curr);
        Side s = Side::north;
        do{
			if (d.hasNbr(s) && visited.count(d.nbr(s)) == 0) {
				// a new edge that we have not assigned an index to

				// fine case
				if (d.hasFineNbr(s)) {
					DomainSignature &nbr_left  = domains.at(d.nbr(s));
					DomainSignature &nbr_right = domains.at(d.nbrRight(s));

					switch (s) {
						case Side::north:
							nbr_left.x_start  = d.x_start;
							nbr_left.y_start  = d.y_start + d.y_length;
							nbr_right.x_start = d.x_start + d.x_length / 2;
							nbr_right.y_start = d.y_start + d.y_length;
							break;
						case Side::east:
							nbr_left.x_start  = d.x_start + d.x_length;
							nbr_left.y_start  = d.y_start + d.y_length / 2;
							nbr_right.x_start = d.x_start + d.x_length;
							nbr_right.y_start = d.y_start;
							break;
						case Side::south:
							nbr_left.x_start  = d.x_start + d.x_length / 2;
							nbr_left.y_start  = d.y_start - d.y_length / 2;
							nbr_right.x_start = d.x_start;
							nbr_right.y_start = d.y_start - d.y_length / 2;
							break;
						case Side::west:
							nbr_left.x_start  = d.x_start - d.x_length / 2;
							nbr_left.y_start  = d.y_start;
							nbr_right.x_start = d.x_start - d.x_length / 2;
							nbr_right.y_start = d.y_start + d.y_length / 2;
					}
					nbr_left.x_length  = d.x_length / 2;
					nbr_left.y_length  = d.y_length / 2;
					nbr_right.x_length = d.x_length / 2;
					nbr_right.y_length = d.y_length / 2;

					// enqueue domains
					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
					if (enqueued.count(d.nbrRight(s)) == 0) {
						queue.push_back(d.nbrRight(s));
						enqueued.insert(d.nbrRight(s));
					}
					// coarse case
				} else if (d.hasCoarseNbr(s)) {
					DomainSignature &nbr = domains.at(d.nbr(s));

					switch (s) {
						case Side::north:
							if (d.leftOfCoarse(s)) {
								nbr.x_start = d.x_start - d.x_length;
								nbr.y_start = d.y_start + d.y_length;
							} else {
								nbr.x_start = d.x_start;
								nbr.y_start = d.y_start + d.y_length;
							}
							break;
                        case Side::east:
							if (d.leftOfCoarse(s)) {
								nbr.x_start = d.x_start + d.x_length;
								nbr.y_start = d.y_start;
							} else {
								nbr.x_start = d.x_start + d.x_length;
								nbr.y_start = d.y_start - d.y_length;
							}
							break;
						case Side::south:
							if (d.leftOfCoarse(s)) {
								nbr.x_start = d.x_start;
								nbr.y_start = d.y_start - d.y_length * 2;
							} else {
								nbr.x_start = d.x_start - d.x_length;
								nbr.y_start = d.y_start - d.y_length * 2;
							}
							break;
						case Side::west:
							if (d.leftOfCoarse(s)) {
								nbr.x_start = d.x_start - d.x_length * 2;
								nbr.y_start = d.y_start - d.y_length;
							} else {
								nbr.x_start = d.x_start - d.x_length * 2;
								nbr.y_start = d.y_start;
							}
					}
					nbr.x_length = 2 * d.x_length;
					nbr.y_length = 2 * d.y_length;

					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
					// normal case
				} else {
					DomainSignature &nbr = domains.at(d.nbr(s));

					switch (s) {
						case Side::north:
							nbr.x_start = d.x_start;
							nbr.y_start = d.y_start + d.y_length;
							break;
						case Side::east:
							nbr.x_start = d.x_start + d.x_length;
							nbr.y_start = d.y_start;
							break;
						case Side::south:
							nbr.x_start = d.x_start;
							nbr.y_start = d.y_start - d.y_length;
							break;
						case Side::west:
							nbr.x_start = d.x_start - d.x_length;
							nbr.y_start = d.y_start;
					}
					nbr.x_length = d.x_length;
					nbr.y_length = d.y_length;

					// enqueue domain
					if (enqueued.count(d.nbr(s)) == 0) {
						queue.push_back(d.nbr(s));
						enqueued.insert(d.nbr(s));
					}
				}
			}
			s++;
		} while (s != Side::north);
	}
    double x_shift = -x_min;
    double y_shift = -y_min;
	double x_scale = x_max - x_min;
	double y_scale = y_max - y_min;
	double scale = x_scale;
    if(y_scale>scale){
        scale=y_scale;
    }
	for (auto &p : domains) {
		DomainSignature &d = p.second;
		d.x_start += x_shift;
		d.y_start += y_shift;
		d.x_start /= scale;
		d.y_start /= scale;
		d.x_length /= scale;
		d.y_length /= scale;
	}
}
DomainSignatureCollection::DomainSignatureCollection(int d_x, int d_y, int rank)
{
	this->rank            = rank;
	num_global_domains    = d_x * d_y;
	num_global_interfaces = 0;
	if (rank == 0) {
		for (int domain_y = 0; domain_y < d_y; domain_y++) {
			for (int domain_x = 0; domain_x < d_x; domain_x++) {
				DomainSignature ds;
				ds.id           = domain_y * d_x + domain_x;
				ds.refine_level = 1;
				if (domain_y != d_y - 1) {
					ds.nbr(Side::north) = (domain_y + 1) * d_x + domain_x;
					ds.proc[0]          = 0;
				}
				if (domain_x != d_x - 1) {
					ds.nbr(Side::east) = domain_y * d_x + domain_x + 1;
					ds.proc[2]         = 0;
				}
				if (domain_y != 0) {
					ds.nbr(Side::south) = (domain_y - 1) * d_x + domain_x;
					ds.proc[4]          = 0;
				}
				if (domain_x != 0) {
					ds.nbr(Side::west) = domain_y * d_x + domain_x - 1;
					ds.proc[6]         = 0;
				}
				ds.x_length    = 1.0 / d_x;
				ds.y_length    = 1.0 / d_y;
				ds.x_start     = 1.0 * domain_x / d_x;
				ds.y_start     = 1.0 * domain_y / d_y;
				domains[ds.id] = ds;
			}
		}
	}
}

DomainSignatureCollection::DomainSignatureCollection(int d_x, int d_y, int rank, bool amr)
{
	this->rank            = rank;
	num_global_domains    = d_x * d_y*5;
	num_global_interfaces = 0;
	if (rank == 0) {
		for (int domain_y = 0; domain_y < d_y; domain_y++) {
			for (int domain_x = 0; domain_x < d_x; domain_x++) {
				DomainSignature ds;
				ds.id = domain_y * d_x + domain_x;
                ds.refine_level=1;
				if (domain_y != d_y - 1) {
					ds.nbr(Side::north) = (domain_y + 1) * d_x + domain_x;
					ds.proc[0]          = 0;
				}
				if (domain_x != d_x - 1) {
					ds.nbr(Side::east) = domain_y * d_x + domain_x + 1;
					ds.proc[2]         = 0;
				}
				if (domain_y != 0) {
					ds.nbr(Side::south) = (domain_y - 1) * d_x + domain_x;
					ds.proc[4]          = 0;
				}
				if (domain_x != 0) {
					ds.nbr(Side::west) = domain_y * d_x + domain_x - 1;
					ds.proc[6] = 0;
				}
				ds.x_length    = 1.0 / d_x;
				ds.y_length    = 1.0 / d_y;
				ds.x_start     = 1.0 * domain_x / d_x;
				ds.y_start     = 1.0 * domain_y / d_y;
				domains[ds.id] = ds;
			}
		}
		// create refined grid
		for (int domain_y = 0; domain_y < d_y * 2; domain_y++) {
			for (int domain_x = 0; domain_x < d_x * 2; domain_x++) {
				DomainSignature ds;
				ds.id           = domain_y * d_x * 2 + domain_x + d_x * d_y;
				ds.refine_level = 2;
				if (domain_y != 2 * d_y - 1) {
					ds.nbr(Side::north) = d_x * d_y + 2 * (domain_y + 1) * d_x + domain_x;
					ds.proc[0]          = 0;
				}
				if (domain_x != 2 * d_x - 1) {
					ds.nbr(Side::east) = d_x * d_y + 2 * domain_y * d_x + domain_x + 1;
					ds.proc[2]         = 0;
				}
				if (domain_y != 0) {
					ds.nbr(Side::south) = d_x * d_y + 2 * (domain_y - 1) * d_x + domain_x;
					ds.proc[4]          = 0;
				}
				if (domain_x != 0) {
					ds.nbr(Side::west) = d_x * d_y + 2 * domain_y * d_x + domain_x - 1;
					ds.proc[6]         = 0;
				}
				ds.x_length    = 1.0 / (2 * d_x);
				ds.y_length    = 1.0 / (2 * d_y);
				ds.x_start     = 1.0 + 1.0 * domain_x / (2 * d_x);
				ds.y_start     = 1.0 * domain_y / (2 * d_y);
				domains[ds.id] = ds;
			}
		}
        //stitch together grids
        for(int i=0;i<d_y;i++){
			DomainSignature &left      = domains[i * d_x + d_x - 1];
			DomainSignature &low_nbr   = domains[d_y * d_x + 2 * i * d_x * 2];
			DomainSignature &high_nbr  = domains[d_y * d_x + (2 * i + 1) * d_x * 2];
			left.nbr(Side::east)       = high_nbr.id;
			left.nbrRight(Side::east)  = low_nbr.id;
			left.nbr_fine[1]           = true;
			low_nbr.nbr(Side::west)    = left.id;
			low_nbr.nbr_coarse[3]      = true;
			low_nbr.left_of_coarse[3]  = false;
			high_nbr.nbr(Side::west)   = left.id;
			high_nbr.nbr_coarse[3]     = true;
			high_nbr.left_of_coarse[3] = true;
		}
	}
}
