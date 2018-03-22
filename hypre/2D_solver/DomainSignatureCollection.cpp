#include "DomainSignatureCollection.h"
#include <iostream>
#include <deque>
#include <set>
#include <fstream>
using namespace std;
DomainSignatureCollection::DomainSignatureCollection(string file_name, int rank)
{
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
	indexInterfacesBFS();
	MPI_Bcast(&num_global_domains, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
void DomainSignatureCollection::divide()
{
	if (rank == 0) {
		std::map<int, DomainSignature> new_domains;
		for (auto p : domains) {
			DomainSignature &z = p.second;
			DomainSignature  a, b, c, d;
			a.id = z.id * 4;
			b.id = z.id * 4 + 1;
			c.id = z.id * 4 + 2;
			d.id = z.id * 4 + 3;

			a.nbr(Side::north) = c.id;
			a.nbr(Side::east)  = b.id;
			b.nbr(Side::north) = d.id;
			b.nbr(Side::west)  = a.id;
			c.nbr(Side::east)  = d.id;
			c.nbr(Side::south) = a.id;
			d.nbr(Side::south) = b.id;
			d.nbr(Side::west)  = c.id;

			// north
			Side s = Side::north;
			if (z.hasFineNbr(s)) {
				c.setHasFineNbr(s);
				c.nbr(s)      = 4 * z.nbr(s);
				c.nbrRight(s) = 4 * z.nbr(s) + 1;

				d.setHasFineNbr(s);
				d.nbr(s)      = 4 * z.nbrRight(s);
				d.nbrRight(s) = 4 * z.nbrRight(s) + 1;
			} else if (z.hasCoarseNbr(s)) {
				c.setHasCoarseNbr(s);
				d.setHasCoarseNbr(s);
				if (z.leftOfCoarse(s)) {
					c.nbr(s) = 4 * z.nbr(s) + 1;
					d.nbr(s) = c.nbr(s);
				} else {
					c.nbr(s) = 4 * z.nbr(s);
					d.nbr(s) = c.nbr(s);
				}
			} else if (z.hasNbr(s)) {
				c.nbr(s) = 4 * z.nbr(s);
				d.nbr(s) = 4 * z.nbr(s) + 1;
			}

			// east
			s = Side::east;
			if (z.hasFineNbr(s)) {
				d.setHasFineNbr(s);
				d.nbr(s)      = 4 * z.nbr(s) + 2;
				d.nbrRight(s) = 4 * z.nbr(s);

				b.setHasFineNbr(s);
				b.nbr(s)      = 4 * z.nbrRight(s) + 2;
				b.nbrRight(s) = 4 * z.nbrRight(s);
			} else if (z.hasCoarseNbr(s)) {
				d.setHasCoarseNbr(s);
				b.setHasCoarseNbr(s);
				if (z.leftOfCoarse(s)) {
					d.nbr(s) = 4 * z.nbr(s);
					b.nbr(s) = d.nbr(s);
				} else {
					d.nbr(s) = 4 * z.nbr(s) + 2;
					b.nbr(s) = d.nbr(s);
				}
			} else if (z.hasNbr(s)) {
				d.nbr(s) = 4 * z.nbr(s) + 2;
				b.nbr(s) = 4 * z.nbr(s);
			}

			// south
			s = Side::south;
			if (z.hasFineNbr(s)) {
				b.setHasFineNbr(s);
				b.nbr(s)      = 4 * z.nbr(s) + 3;
				b.nbrRight(s) = 4 * z.nbr(s) + 2;

				a.setHasFineNbr(s);
				a.nbr(s)      = 4 * z.nbrRight(s) + 3;
				a.nbrRight(s) = 4 * z.nbrRight(s) + 2;
			} else if (z.hasCoarseNbr(s)) {
				b.setHasCoarseNbr(s);
				a.setHasCoarseNbr(s);
				if (z.leftOfCoarse(s)) {
					b.nbr(s) = 4 * z.nbr(s) + 2;
					a.nbr(s) = b.nbr(s);
				} else {
					b.nbr(s) = 4 * z.nbr(s) + 3;
					a.nbr(s) = b.nbr(s);
				}
			} else if (z.hasNbr(s)) {
				b.nbr(s) = 4 * z.nbr(s) + 3;
				a.nbr(s) = 4 * z.nbr(s) + 2;
			}

			// west
			s = Side::west;
			if (z.hasFineNbr(s)) {
				a.setHasFineNbr(s);
				a.nbr(s)      = 4 * z.nbr(s) + 1;
				a.nbrRight(s) = 4 * z.nbr(s) + 3;

				c.setHasFineNbr(s);
				c.nbr(s)      = 4 * z.nbrRight(s) + 1;
				c.nbrRight(s) = 4 * z.nbrRight(s) + 3;
			} else if (z.hasCoarseNbr(s)) {
				a.setHasCoarseNbr(s);
				c.setHasCoarseNbr(s);
				if (z.leftOfCoarse(s)) {
					a.nbr(s) = 4 * z.nbr(s) + 3;
					c.nbr(s) = a.nbr(s);
				} else {
					a.nbr(s) = 4 * z.nbr(s) + 1;
					c.nbr(s) = a.nbr(s);
				}
			} else if (z.hasNbr(s)) {
				a.nbr(s) = 4 * z.nbr(s) + 1;
				c.nbr(s) = 4 * z.nbr(s) + 3;
			}

			// add them
			new_domains[a.id] = a;
			new_domains[b.id] = b;
			new_domains[c.id] = c;
			new_domains[d.id] = d;
		}

		// replace old domains
		domains = new_domains;

		determineCoarseness();
		determineAmrLevel();
		determineXY();
		num_global_domains = domains.size();
	}
	indexInterfacesBFS();
	MPI_Bcast(&num_global_domains, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
void DomainSignatureCollection::determineCoarseness()
{
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
		do {
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
		indexInterfacesBFS();
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
		indexInterfacesBFS();
}
void DomainSignatureCollection::indexInterfacesBFS()
{
	int curr_i = 0;
    if(domains.size()!=0){
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
			if (d.hasNbr(s) && d.index(s) == -1) {
				// a new edge that we have not assigned an index to
				d.index(s) = curr_i;
				curr_i++;

				// fine case
				if (d.hasFineNbr(s)) {
					DomainSignature &nbr_left  = domains.at(d.nbr(s));
					DomainSignature &nbr_right = domains.at(d.nbrRight(s));

					// set center indexes
					nbr_left.indexCenter(!s)  = d.index(s);
					nbr_right.indexCenter(!s) = d.index(s);

					// set left and right indexes index
					nbr_left.index(!s) = curr_i;
					curr_i++;
					nbr_right.index(!s) = curr_i;
					curr_i++;

                    // set refined indexes
					d.indexRefinedLeft(s)  = nbr_left.index(!s);
					d.indexRefinedRight(s) = nbr_right.index(!s);

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
					DomainSignature &nbr      = domains.at(d.nbr(s));
					int              buddy_id = -1;
					if (d.leftOfCoarse(s)) {
						DomainSignature &buddy = domains.at(nbr.nbrRight(!s));
						buddy_id               = buddy.id;

						nbr.indexRefinedLeft(!s) = d.index(s);

						nbr.indexRefinedRight(!s) = curr_i;
						buddy.index(s)            = curr_i;
						curr_i++;

						d.indexCenter(s)     = curr_i;
						nbr.index(!s)        = curr_i;
						buddy.indexCenter(s) = curr_i;
						curr_i++;
					} else {
						DomainSignature &buddy = domains.at(nbr.nbr(!s));
						buddy_id               = buddy.id;

						nbr.indexRefinedRight(!s) = d.index(s);

						nbr.indexRefinedLeft(!s) = curr_i;
						buddy.index(s)           = curr_i;
						curr_i++;

						d.indexCenter(s)     = curr_i;
						nbr.index(!s)        = curr_i;
						buddy.indexCenter(s) = curr_i;
						curr_i++;
					}

					// enqueue domains
					if (enqueued.count(nbr.id) == 0) {
						queue.push_back(nbr.id);
						enqueued.insert(nbr.id);
					}
					if (enqueued.count(buddy_id) == 0) {
						queue.push_back(buddy_id);
						enqueued.insert(buddy_id);
					}
					// normal case
				} else {
					DomainSignature &nbr = domains.at(d.nbr(s));
					nbr.index(!s)        = d.index(s);
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
	num_global_interfaces = curr_i;
	matrix_j_low          = 0;
	matrix_j_high         = num_global_interfaces;
	MPI_Bcast(&num_global_interfaces, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
void DomainSignatureCollection::zoltanBalance()
{
	int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
	matrix_j_low = num_global_interfaces * rank / size;
	matrix_j_high = num_global_interfaces * (rank + 1) / size;
    cerr << matrix_j_low << "," << matrix_j_high << "," << rank<<","<<size<< endl;
	Zoltan *zz = new Zoltan(MPI_COMM_WORLD);

	// parameters
	zz->Set_Param("LB_METHOD", "GRAPH");       /* Zoltan method: "BLOCK" */
	zz->Set_Param("LB_APPROACH", "PARTITION"); /* Zoltan method: "BLOCK" */
	zz->Set_Param("NUM_GID_ENTRIES", "1");     /* global ID is 1 integer */
	zz->Set_Param("NUM_LID_ENTRIES", "1");     /* local ID is 1 integer */
	zz->Set_Param("OBJ_WEIGHT_DIM", "0");      /* we omit object weights */
	zz->Set_Param("AUTO_MIGRATE", "TRUE");     /* we omit object weights */

	// Query functions
	zz->Set_Num_Obj_Fn(DomainSignatureCollection::get_number_of_objects, this);
	zz->Set_Obj_List_Fn(DomainSignatureCollection::get_object_list, this);
	zz->Set_Pack_Obj_Multi_Fn(DomainSignatureCollection::pack_objects, this);
	zz->Set_Unpack_Obj_Multi_Fn(DomainSignatureCollection::unpack_objects, this);
	zz->Set_Obj_Size_Multi_Fn(DomainSignatureCollection::object_sizes, this);
	zz->Set_Num_Edges_Fn(DomainSignatureCollection::numInterfaces, this);
	zz->Set_Edge_List_Fn(DomainSignatureCollection::interfaceList, this);
	// zz->Set_Geom_Fn(DomainSignatureCollection::coord, this);
	// zz->Set_Num_Geom_Fn(DomainSignatureCollection::dimensions, this);

	////////////////////////////////////////////////////////////////
	// Zoltan can now partition the objects in this collection.
	// In this simple example, we assume the number of partitions is
	// equal to the number of processes.  Process rank 0 will own
	// partition 0, process rank 1 will own partition 1, and so on.
	////////////////////////////////////////////////////////////////

	int           changes;
	int           numGidEntries;
	int           numLidEntries;
	int           numImport;
	ZOLTAN_ID_PTR importGlobalIds;
	ZOLTAN_ID_PTR importLocalIds;
	int *         importProcs;
	int *         importToPart;
	int           numExport;
	ZOLTAN_ID_PTR exportGlobalIds;
	ZOLTAN_ID_PTR exportLocalIds;
	int *         exportProcs;
	int *         exportToPart;

	int rc = zz->LB_Partition(changes, numGidEntries, numLidEntries, numImport, importGlobalIds,
	                          importLocalIds, importProcs, importToPart, numExport, exportGlobalIds,
	                          exportLocalIds, exportProcs, exportToPart);

	if (rc != ZOLTAN_OK) {
		cerr << "zoltan error\n";
		delete zz;
		exit(0);
	}
	cout << "I have " << domains.size() << " domains: ";

    int prev=-100;
    bool range = false;
	for (auto &p : domains)
	{
        int curr = p.second.id;
		if (curr != prev + 1 && !range) {
			cout << curr << "-";
            range = true;
		}else if(curr != prev + 1 && range){
            cout << prev << " " << curr << "-";
        }
		prev = curr;
        
	}

	cout <<prev<< "\n";
}

/*int DomainSignatureCollection::dimensions(void *data, int *ierr)
{
	*ierr = ZOLTAN_OK;
	return 2;
}*/
/*void DomainSignatureCollection::coord(void *data, int num_gid_entries, int num_lid_entries,
                                      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                                      double *geom_vec, int *ierr){
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;

	auto &ds    = dsc->domains[*global_id];
	geom_vec[0] = ds.id % dsc->d_y;
	geom_vec[1] = ds.id % dsc->d_x;
}*/
// query functions that respond to requests from Zoltan
int DomainSignatureCollection::get_number_of_objects(void *data, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;

	return dsc->domains.size();
}

void DomainSignatureCollection::get_object_list(void *data, int sizeGID, int sizeLID,
                                                ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                                int wgt_dim, float *obj_wgts, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;

	// In this example, return the IDs of our objects, but no weights.
	// Zoltan will assume equally weighted objects.

	int i = 0;
	for (auto &p : dsc->domains) {
		globalID[i] = p.first;
		localID[i]  = p.first;
		i++;
	}
	return;
}

void DomainSignatureCollection::object_sizes(void *data, int num_gid_entries, int num_lid_entries,
                                             int num_ids, ZOLTAN_ID_PTR global_ids,
                                             ZOLTAN_ID_PTR local_ids, int *sizes, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	for (int i = 0; i < num_ids; i++) {
		sizes[i] = sizeof(dsc->domains[global_ids[i]]);
	}
}
void DomainSignatureCollection::pack_objects(void *data, int num_gid_entries, int num_lid_entries,
                                             int num_ids, ZOLTAN_ID_PTR global_ids,
                                             ZOLTAN_ID_PTR local_ids, int *dest, int *sizes,
                                             int *idx, char *buf, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	for (int i = 0; i < num_ids; i++) {
		auto &ds = dsc->domains[global_ids[i]];
		for (int q = 0; q < 8; q++) {
			if (ds.nbr_id[q] != -1) {
				auto &nbr   = dsc->domains[ds.nbr_id[q]];
				nbr.proc[q] = dest[i];
			}
		}
	}
	for (int i = 0; i < num_ids; i++) {
		*((DomainSignature *) &buf[idx[i]]) = dsc->domains[global_ids[i]];
		dsc->domains.erase(global_ids[i]);
	}
}
void DomainSignatureCollection::unpack_objects(void *data, int num_gid_entries, int num_ids,
                                               ZOLTAN_ID_PTR global_ids, int *sizes, int *idx,
                                               char *buf, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	for (int i = 0; i < num_ids; i++) {
		dsc->domains[global_ids[i]] = *((DomainSignature *) &buf[idx[i]]);
	}
}
int DomainSignatureCollection::numInterfaces(void *data, int num_gid_entries, int num_lid_entries,
                                             ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                                             int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	auto &ds                       = dsc->domains[*global_id];
	int   num_iface                = 0;
	for (int q = 0; q < 8; q++) {
		if (ds.nbr_id[q] != -1) num_iface++;
	}
	return num_iface;
}
void DomainSignatureCollection::interfaceList(void *data, int num_gid_entries, int num_lid_entries,
                                              ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                                              ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs,
                                              int wgt_dim, float *ewgts, int *ierr)
{
	DomainSignatureCollection *dsc = (DomainSignatureCollection *) data;
	*ierr                          = ZOLTAN_OK;
	auto &ds                       = dsc->domains[*global_id];
	int   i                        = 0;
	for (int q = 0; q < 8; q++) {
		if (ds.nbr_id[q] != -1) {
			nbor_global_id[i] = ds.nbr_id[q];
			nbor_procs[i]     = ds.proc[q];
			i++;
		}
	}
}
