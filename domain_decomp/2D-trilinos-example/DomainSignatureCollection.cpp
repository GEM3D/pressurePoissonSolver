#include "DomainSignatureCollection.h"
#include <deque>
#include <fstream>
#include <iostream>
#include <set>
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
			if (ds.id != -1) {
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
		}
		determineCoarseness();
		determineAmrLevel();
		determineXY();
		num_global_domains = domains.size();
	}
	indexInterfacesBFS();
	MPI_Bcast(&num_global_domains, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
DomainSignatureCollection::DomainSignatureCollection(int d_x, int d_y, int rank)
{
	this->rank         = rank;
	num_global_domains = d_x * d_y;
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
void DomainSignatureCollection::enumerateIfaces()
{
	map<int, Iface> ifaces_new;
	ifaces = ifaces_new;
	for (auto p : domains) {
		DomainSignature &ds = p.second;
		Side             s  = Side::north;
		do {
			int index = ds.globalIndex(s);
			if (index != -1 && ifaces.count(index) == 0) {
				Iface &iface = ifaces[index];
				iface.id     = index;
				iface.y_axis = s == Side::west || s == Side::east;
				bool left    = s == Side::north || s == Side::east;
				if (ds.hasFineNbr(s)) {
					if (left) {
						iface.type  = IfaceType::coarse_on_left;
						iface.left  = ds;
						iface.right = domains[ds.nbr(s)];
						iface.extra = domains[ds.nbrRight(s)];
					} else {
						iface.type  = IfaceType::coarse_on_right;
						iface.right = ds;
						iface.left  = domains[ds.nbr(s)];
						iface.extra = domains[ds.nbrRight(s)];
					}
				} else if (ds.hasCoarseNbr(s)) {
					if (left) {
						if (ds.leftOfCoarse(s)) {
							iface.type = IfaceType::refined_on_left_left_of_coarse;
						} else {
							iface.type = IfaceType::refined_on_left_right_of_coarse;
						}
						iface.left  = ds;
						iface.right = domains[ds.nbr(s)];
					} else {
						if (ds.leftOfCoarse(s)) {
							iface.type = IfaceType::refined_on_right_left_of_coarse;
						} else {
							iface.type = IfaceType::refined_on_right_right_of_coarse;
						}
						iface.right = ds;
						iface.left  = domains[ds.nbr(s)];
					}
				} else {
					iface.type = IfaceType::normal;
					if (left) {
						iface.left  = ds;
						iface.right = domains[ds.nbr(s)];
					} else {
						iface.right = ds;
						iface.left  = domains[ds.nbr(s)];
					}
				}
			}
			s++;
		} while (s != Side::north);
	}
	num_pins = 0;
	for (auto p : ifaces) {
		num_pins += p.second.getPins().size();
	}
	indexIfacesLocal();
	indexDomainIfacesLocal();
}
DomainSignatureCollection::DomainSignatureCollection(int d_x, int d_y, int rank, bool amr)
{
	this->rank         = rank;
	num_global_domains = d_x * d_y * 5;
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
		// stitch together grids
		for (int i = 0; i < d_y; i++) {
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
	enumerateIfaces();
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
void DomainSignatureCollection::determineAmrLevel()
{
	set<int>   visited;
	set<int>   enqueued;
	deque<int> queue;
	int        first = domains.begin()->first;
	queue.push_back(first);
	enqueued.insert(first);
	int min_level = 1;
	while (!queue.empty()) {
		int              curr       = queue.front();
		DomainSignature &d          = domains.at(curr);
		int              curr_level = d.refine_level;
		queue.pop_front();
		visited.insert(curr);
		Side s = Side::north;
		do {
			if (d.hasNbr(s) && visited.count(d.nbr(s)) == 0) {
				// fine case
				if (d.hasFineNbr(s)) {
					DomainSignature &nbr_left  = domains.at(d.nbr(s));
					DomainSignature &nbr_right = domains.at(d.nbrRight(s));

					nbr_left.refine_level  = curr_level + 1;
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

					nbr.refine_level = curr_level - 1;
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
void DomainSignatureCollection::determineXY()
{
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
		if (d.x_start + d.x_length > x_max) {
			x_max = d.x_start + d.x_length;
		}
		if (d.y_start + d.y_length > y_max) {
			y_max = d.y_start + d.y_length;
		}
		visited.insert(curr);
		Side s = Side::north;
		do {
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
	double scale   = x_scale;
	if (y_scale > scale) {
		scale = y_scale;
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
void DomainSignatureCollection::indexInterfacesBFS()
{
	int curr_i = 0;
	if (domains.size() != 0) {
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
				if (d.hasNbr(s) && d.globalIndex(s) == -1) {
					// a new edge that we have not assigned an index to
					d.globalIndex(s) = curr_i;
					curr_i++;

					// fine case
					if (d.hasFineNbr(s)) {
						DomainSignature &nbr_left  = domains.at(d.nbr(s));
						DomainSignature &nbr_right = domains.at(d.nbrRight(s));

						// set center indexes
						nbr_left.globalIndexCenter(!s)  = d.globalIndex(s);
						nbr_right.globalIndexCenter(!s) = d.globalIndex(s);

						// set left and right indexes index
						nbr_left.globalIndex(!s) = curr_i;
						curr_i++;
						nbr_right.globalIndex(!s) = curr_i;
						curr_i++;

						// set refined indexes
						d.globalIndexRefinedLeft(s)  = nbr_left.globalIndex(!s);
						d.globalIndexRefinedRight(s) = nbr_right.globalIndex(!s);

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

							nbr.globalIndexRefinedLeft(!s) = d.globalIndex(s);

							nbr.globalIndexRefinedRight(!s) = curr_i;
							buddy.globalIndex(s)            = curr_i;
							curr_i++;

							d.globalIndexCenter(s)     = curr_i;
							nbr.globalIndex(!s)        = curr_i;
							buddy.globalIndexCenter(s) = curr_i;
							curr_i++;
						} else {
							DomainSignature &buddy = domains.at(nbr.nbr(!s));
							buddy_id               = buddy.id;

							nbr.globalIndexRefinedRight(!s) = d.globalIndex(s);

							nbr.globalIndexRefinedLeft(!s) = curr_i;
							buddy.globalIndex(s)           = curr_i;
							curr_i++;

							d.globalIndexCenter(s)     = curr_i;
							nbr.globalIndex(!s)        = curr_i;
							buddy.globalIndexCenter(s) = curr_i;
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
						nbr.globalIndex(!s)  = d.globalIndex(s);
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
	enumerateIfaces();
}
void DomainSignatureCollection::zoltanBalance()
{
	zoltanBalanceIfaces();
	zoltanBalanceDomains();
	indexIfacesLocal();
	indexDomainIfacesLocal();
}

void DomainSignatureCollection::zoltanBalanceIfaces()
{
	Zoltan *zz = new Zoltan(MPI_COMM_WORLD);

	// parameters
	zz->Set_Param("LB_METHOD", "HYPERGRAPH");
	zz->Set_Param("LB_APPROACH", "PARTITION");
	zz->Set_Param("NUM_GID_ENTRIES", "1");
	zz->Set_Param("NUM_LID_ENTRIES", "1");
	zz->Set_Param("OBJ_WEIGHT_DIM", "1");
	zz->Set_Param("AUTO_MIGRATE", "TRUE");
	zz->Set_Param("PHG_EDGE_SIZE_THRESHOLD", "1.0");
	zz->Set_Param("IMBALANCE_TOL[0]", "1.0");

	// Query functions
	zz->Set_Num_Obj_Fn(IfaceZoltanHelper::get_number_of_objects, this);
	zz->Set_Obj_List_Fn(IfaceZoltanHelper::get_object_list, this);
	zz->Set_Pack_Obj_Multi_Fn(IfaceZoltanHelper::pack_objects, this);
	zz->Set_Unpack_Obj_Multi_Fn(IfaceZoltanHelper::unpack_objects, this);
	zz->Set_Obj_Size_Multi_Fn(IfaceZoltanHelper::object_sizes, this);
	zz->Set_HG_Size_CS_Fn(IfaceZoltanHelper::ZOLTAN_HG_SIZE_CS_FN, this);
	zz->Set_HG_CS_Fn(IfaceZoltanHelper::ZOLTAN_HG_CS_FN, this);
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
	delete zz;
	cerr << "I have " << ifaces.size() << " ifaces: ";

#if DD_DEBUG
	int  prev  = -100;
	bool range = false;
	for (auto &p : ifaces) {
		int curr = p.second.id;
		if (curr != prev + 1 && !range) {
			cout << curr << "-";
			range = true;
		} else if (curr != prev + 1 && range) {
			cout << prev << " " << curr << "-";
		}
		prev = curr;
	}

	cout << prev << "\n";
#endif
	cerr << endl;
}
void DomainSignatureCollection::zoltanBalanceDomains()
{
	Zoltan *zz = new Zoltan(MPI_COMM_WORLD);

	// parameters
	zz->Set_Param("LB_METHOD", "GRAPH");       /* Zoltan method: "BLOCK" */
	zz->Set_Param("LB_APPROACH", "PARTITION"); /* Zoltan method: "BLOCK" */
	zz->Set_Param("NUM_GID_ENTRIES", "1");     /* global ID is 1 integer */
	zz->Set_Param("NUM_LID_ENTRIES", "1");     /* local ID is 1 integer */
	zz->Set_Param("OBJ_WEIGHT_DIM", "0");      /* we omit object weights */
	zz->Set_Param("AUTO_MIGRATE", "TRUE");     /* we omit object weights */

	// Query functions
	zz->Set_Num_Obj_Fn(DomainZoltanHelper::get_number_of_objects, this);
	zz->Set_Obj_List_Fn(DomainZoltanHelper::get_object_list, this);
	zz->Set_Pack_Obj_Multi_Fn(DomainZoltanHelper::pack_objects, this);
	zz->Set_Unpack_Obj_Multi_Fn(DomainZoltanHelper::unpack_objects, this);
	zz->Set_Obj_Size_Multi_Fn(DomainZoltanHelper::object_sizes, this);
	zz->Set_Num_Edges_Fn(DomainZoltanHelper::numInterfaces, this);
	zz->Set_Edge_List_Fn(DomainZoltanHelper::interfaceList, this);
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

#if DD_DEBUG
	int  prev  = -100;
	bool range = false;
	for (auto &p : domains) {
		int curr = p.second.id;
		if (curr != prev + 1 && !range) {
			cout << curr << "-";
			range = true;
		} else if (curr != prev + 1 && range) {
			cout << prev << " " << curr << "-";
		}
		prev = curr;
	}

	cout << prev << "\n";
#endif
	cout << endl;
}

void DomainSignatureCollection::indexIfacesLocal()
{
	vector<int> map_vec;
	vector<int> off_proc_map_vec;
	map<int, int> rev_map;
	if (!ifaces.empty()) {
		deque<int> queue;
		set<int>   enqueued;
		queue.push_back(ifaces.begin()->first);
		enqueued.insert(ifaces.begin()->first);
		deque<int> off_proc_ifaces;
		int        curr_i = 0;
		while (!queue.empty()) {
			int i = queue.front();
			queue.pop_front();
			map_vec.push_back(i);
			Iface &iface       = ifaces[i];
			rev_map[i]         = curr_i;
			ifaces[i].id_local = curr_i;
			curr_i++;
			set<MatrixBlock> blocks = iface.getGlobalRowBlocks();
			for (MatrixBlock b : blocks) {
				if (!enqueued.count(b.j)) {
					enqueued.insert(b.j);
					if (ifaces.count(b.j)) {
						queue.push_back(b.j);
					} else {
						off_proc_map_vec.push_back(b.j);
					}
				}
			}
		}
		// map off proc
		for (int i : off_proc_map_vec) {
			rev_map[i] = curr_i;
			curr_i++;
		}
		for (auto &p : ifaces) {
			p.second.setLocalIndexes(rev_map);
		}
	}
	iface_rev_map          = rev_map;
	iface_map_vec          = map_vec;
	iface_off_proc_map_vec = off_proc_map_vec;
}
void DomainSignatureCollection::indexDomainIfacesLocal()
{
	vector<int> map_vec;
	map<int, int> rev_map;
	if (!domains.empty()) {
		deque<int> queue;
		set<int>   enqueued;
		queue.push_back(domains.begin()->first);
		enqueued.insert(domains.begin()->first);
		deque<int> off_proc_ifaces;
		int        curr_i = 0;
		while (!queue.empty()) {
			int i = queue.front();
			queue.pop_front();
			DomainSignature &ds = domains[i];
			for (int global_i : ds.global_i) {
				if (global_i != -1 && rev_map.count(global_i) == 0) {
					rev_map[global_i] = curr_i;
					map_vec.push_back(global_i);
					curr_i++;
				}
			}
			for (int global_i : ds.global_i_refined) {
				if (global_i != -1 && rev_map.count(global_i) == 0) {
					rev_map[global_i] = curr_i;
					map_vec.push_back(global_i);
					curr_i++;
				}
			}
			for (int global_i : ds.global_i_center) {
				if (global_i != -1 && rev_map.count(global_i) == 0) {
					rev_map[global_i] = curr_i;
					map_vec.push_back(global_i);
					curr_i++;
				}
			}
			for (int nbr : ds.nbr_id) {
				if (nbr != -1 && enqueued.count(nbr) == 0) {
					enqueued.insert(nbr);
					if (domains.count(nbr)) {
						queue.push_back(nbr);
					}
				}
			}
		}
		for (auto &p : domains) {
			p.second.setLocalIndexes(rev_map);
		}
	}
	domain_rev_map = rev_map;
	domain_map_vec = map_vec;
}
std::set<MatrixBlock> Iface::getGlobalRowBlocks()
{
	std::set<MatrixBlock> blocks;
	auto getBlockType = [=](bool left) {
		BlockType ret = BlockType::plain;
		switch (type) {
			case IfaceType::normal:
				ret = BlockType::plain;
				break;
			case IfaceType::coarse_on_left:
				if (left) {
					ret = BlockType::coarse;
				} else {
					ret = BlockType::fine_out_left;
				}
				break;
			case IfaceType::coarse_on_right:
				if (left) {
					ret = BlockType::fine_out_left;
				} else {
					ret = BlockType::coarse;
				}
				break;
			case IfaceType::refined_on_left_left_of_coarse:
				if (left) {
					ret = BlockType::fine;
				} else {
					ret = BlockType::coarse_out_left;
				}
				break;
			case IfaceType::refined_on_left_right_of_coarse:
				if (left) {
					ret = BlockType::fine;
				} else {
					ret = BlockType::coarse_out_right;
				}
				break;
			case IfaceType::refined_on_right_left_of_coarse:
				if (left) {
					ret = BlockType::coarse_out_left;
				} else {
					ret = BlockType::fine;
				}
				break;
			case IfaceType::refined_on_right_right_of_coarse:
				if (left) {
					ret = BlockType::coarse_out_right;
				} else {
					ret = BlockType::fine;
				}
				break;
		}
		return ret;
	};
	int i = id;
	// left domain(s)
	{
		Side iface_s = Side::north;
		if (y_axis) {
			iface_s = Side::east;
		}
		// north block
		{
			Side           s    = Side::north;
			Side           main = iface_s + s;
			std::bitset<4> neumann;
			int            j = left.globalIndex(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, false, getBlockType(true));
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = left.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, getBlockType(true));
				blocks.insert(b);
			}
		}
		// south block
		{
			Side s    = Side::south;
			Side main = iface_s + s;
			int  j    = left.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, getBlockType(true));
				blocks.insert(b);
			}
		}
		// west block
		{
			Side s    = Side::east;
			Side rel  = Side::west;
			Side main = iface_s + rel;
			int  j    = left.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, getBlockType(true));
				blocks.insert(b);
			}
		}
	}
	// right domain(s)
	{
		Side iface_s = Side::south;
		if (y_axis) {
			iface_s = Side::west;
		}
		// north block
		{
			Side           s    = Side::north;
			Side           main = iface_s + s;
			std::bitset<4> neumann;
			int            j = right.globalIndex(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, false, getBlockType(false));
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = right.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, getBlockType(false));
				blocks.insert(b);
			}
		}
		// south block
		{
			Side s    = Side::south;
			Side main = iface_s + s;
			int  j    = right.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, getBlockType(false));
				blocks.insert(b);
			}
		}
		// west block
		{
			Side s    = Side::east;
			Side rel  = Side::west;
			Side main = iface_s + rel;
			int  j    = right.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, getBlockType(false));
				blocks.insert(b);
			}
		}
	}
	if (type == IfaceType::coarse_on_left || type == IfaceType::coarse_on_right) {
		Side iface_s = Side::south;
		if (y_axis) {
			iface_s = Side::west;
		}
		if (type == IfaceType::coarse_on_right) {
			iface_s = Side::north;
			if (y_axis) {
				iface_s = Side::east;
			}
		}
		// north block
		{
			Side           s    = Side::north;
			Side           main = iface_s + s;
			std::bitset<4> neumann;
			int            j = extra.globalIndex(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, false, BlockType::fine_out_right);
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = extra.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, BlockType::fine_out_right);
				blocks.insert(b);
			}
		}
		// south block
		{
			Side s    = Side::south;
			Side main = iface_s + s;
			int  j    = extra.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, BlockType::fine_out_right);
				blocks.insert(b);
			}
		}
		// west block
		{
			Side s    = Side::east;
			Side rel  = Side::west;
			Side main = iface_s + rel;
			int  j    = extra.globalIndex(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, BlockType::fine_out_right);
				blocks.insert(b);
			}
		}
	}
	return blocks;
}
std::set<MatrixBlock> Iface::getGlobalColBlocks()
{
	std::set<MatrixBlock> blocks;
	auto addBlocks = [&](DomainSignature &ds, Side main, Side s) {
		std::bitset<4> neumann = ds.neumannRelative(main);
		int            j       = ds.globalIndex(main);
		bool zp = ds.zero_patch;
		if (ds.hasFineNbr(main + s)) {
			MatrixBlock b(ds.globalIndex(main + s), j, main, s, neumann, zp, BlockType::coarse);
			MatrixBlock c(ds.globalIndexRefinedLeft(main + s), j, main, s, neumann, zp,
			              BlockType::coarse_out_left);
			MatrixBlock d(ds.globalIndexRefinedRight(main + s), j, main, s, neumann, zp,
			              BlockType::coarse_out_right);
			blocks.insert(b);
			blocks.insert(c);
			blocks.insert(d);
		} else if (ds.hasCoarseNbr(main + s)) {
			MatrixBlock b(ds.globalIndex(main + s), j, main, s, neumann, zp, BlockType::fine);
			blocks.insert(b);
			if (ds.leftOfCoarse(main + s)) {
				MatrixBlock c(ds.globalIndexCenter(main + s), j, main, s, neumann, zp,
				              BlockType::fine_out_left);
				blocks.insert(c);
			} else {
				MatrixBlock c(ds.globalIndexCenter(main + s), j, main, s, neumann, zp,
				              BlockType::fine_out_right);
				blocks.insert(c);
			}
		} else {
			MatrixBlock b(ds.globalIndex(main + s), j, main, s, neumann, zp, BlockType::plain);
			blocks.insert(b);
		}
	};
	// left block
	{
		Side iface_s = Side::north;
		if (y_axis) {
			iface_s = Side::east;
		}
		if (left.globalIndex(iface_s) == id) {
			// north block
			{
				Side s = Side::north;
				addBlocks(left, iface_s, s);
			}
			// east block
			{
				Side s = Side::east;
				if (left.hasNbr(iface_s + s)) {
					addBlocks(left, iface_s, s);
				}
			}
			// south block
			{
				Side s = Side::south;
				if (left.hasNbr(iface_s + s)) {
					addBlocks(left, iface_s, s);
				}
			}
			// west block
			{
				Side s = Side::west;
				if (left.hasNbr(iface_s + s)) {
					addBlocks(left, iface_s, s);
				}
			}
		}
	}
	// right block
	{
		Side iface_s = Side::south;
		if (y_axis) {
			iface_s = Side::west;
		}
		if (right.globalIndex(iface_s) == id) {
			// north block
			{
				Side s = Side::north;
				addBlocks(right, iface_s, s);
			}
			// east block
			{
				Side s = Side::east;
				if (right.hasNbr(iface_s + s)) {
					addBlocks(right, iface_s, s);
				}
			}
			// south block
			{
				Side s = Side::south;
				if (right.hasNbr(iface_s + s)) {
					addBlocks(right, iface_s, s);
				}
			}
			// west block
			{
				Side s = Side::west;
				if (right.hasNbr(iface_s + s)) {
					addBlocks(right, iface_s, s);
				}
			}
		}
	}

	return blocks;
}
//*******
// local
//*******
std::set<MatrixBlock> Iface::getRowBlocks()
{
	std::set<MatrixBlock> blocks;
	auto getBlockType = [=](bool left) {
		BlockType ret = BlockType::plain;
		switch (type) {
			case IfaceType::normal:
				ret = BlockType::plain;
				break;
			case IfaceType::coarse_on_left:
				if (left) {
					ret = BlockType::coarse;
				} else {
					ret = BlockType::fine_out_left;
				}
				break;
			case IfaceType::coarse_on_right:
				if (left) {
					ret = BlockType::fine_out_left;
				} else {
					ret = BlockType::coarse;
				}
				break;
			case IfaceType::refined_on_left_left_of_coarse:
				if (left) {
					ret = BlockType::fine;
				} else {
					ret = BlockType::coarse_out_left;
				}
				break;
			case IfaceType::refined_on_left_right_of_coarse:
				if (left) {
					ret = BlockType::fine;
				} else {
					ret = BlockType::coarse_out_right;
				}
				break;
			case IfaceType::refined_on_right_left_of_coarse:
				if (left) {
					ret = BlockType::coarse_out_left;
				} else {
					ret = BlockType::fine;
				}
				break;
			case IfaceType::refined_on_right_right_of_coarse:
				if (left) {
					ret = BlockType::coarse_out_right;
				} else {
					ret = BlockType::fine;
				}
				break;
		}
		return ret;
	};
	int i = id_local;
	// left domain(s)
	{
		bool zp      = left.zero_patch;
		Side iface_s = Side::north;
		if (y_axis) {
			iface_s = Side::east;
		}
		// north block
		{
			Side           s       = Side::north;
			Side           main    = iface_s + s;
			std::bitset<4> neumann = left.neumannRelative(main);
			int            j       = left.index(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, zp, getBlockType(true));
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = left.index(main);
			if (j != -1) {
				std::bitset<4> neumann = left.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, getBlockType(true));
				blocks.insert(b);
			}
		}
		// south block
		{
			Side s    = Side::south;
			Side main = iface_s + s;
			int  j    = left.index(main);
			if (j != -1) {
				std::bitset<4> neumann = left.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, getBlockType(true));
				blocks.insert(b);
			}
		}
		// west block
		{
			Side s    = Side::east;
			Side rel  = Side::west;
			Side main = iface_s + rel;
			int  j    = left.index(main);
			if (j != -1) {
				std::bitset<4> neumann = left.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, getBlockType(true));
				blocks.insert(b);
			}
		}
	}
	// right domain(s)
	{
		bool zp      = right.zero_patch;
		Side iface_s = Side::south;
		if (y_axis) {
			iface_s = Side::west;
		}
		// north block
		{
			Side           s       = Side::north;
			Side           main    = iface_s + s;
			std::bitset<4> neumann = right.neumannRelative(main);
			int            j       = right.index(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, zp, getBlockType(false));
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = right.index(main);
			if (j != -1) {
				std::bitset<4> neumann = right.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, getBlockType(false));
				blocks.insert(b);
			}
		}
		// south block
		{
			Side s    = Side::south;
			Side main = iface_s + s;
			int  j    = right.index(main);
			if (j != -1) {
				std::bitset<4> neumann = right.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, getBlockType(false));
				blocks.insert(b);
			}
		}
		// west block
		{
			Side s    = Side::east;
			Side rel  = Side::west;
			Side main = iface_s + rel;
			int  j    = right.index(main);
			if (j != -1) {
				std::bitset<4> neumann = right.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, getBlockType(false));
				blocks.insert(b);
			}
		}
	}
	if (type == IfaceType::coarse_on_left || type == IfaceType::coarse_on_right) {
		bool zp      = extra.zero_patch;
		Side iface_s = Side::south;
		if (y_axis) {
			iface_s = Side::west;
		}
		if (type == IfaceType::coarse_on_right) {
			iface_s = Side::north;
			if (y_axis) {
				iface_s = Side::east;
			}
		}
		// north block
		{
			Side           s       = Side::north;
			Side           main    = iface_s + s;
			std::bitset<4> neumann = extra.neumannRelative(main);
			int            j       = extra.index(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, zp, BlockType::fine_out_right);
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = extra.index(main);
			if (j != -1) {
				std::bitset<4> neumann = extra.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, BlockType::fine_out_right);
				blocks.insert(b);
			}
		}
		// south block
		{
			Side s    = Side::south;
			Side main = iface_s + s;
			int  j    = extra.index(main);
			if (j != -1) {
				std::bitset<4> neumann = extra.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, BlockType::fine_out_right);
				blocks.insert(b);
			}
		}
		// west block
		{
			Side s    = Side::east;
			Side rel  = Side::west;
			Side main = iface_s + rel;
			int  j    = extra.index(main);
			if (j != -1) {
				std::bitset<4> neumann = extra.neumannRelative(main);
				MatrixBlock    b(i, j, main, s, neumann, zp, BlockType::fine_out_right);
				blocks.insert(b);
			}
		}
	}
	return blocks;
}
std::set<MatrixBlock> Iface::getColBlocks()
{
	std::set<MatrixBlock> blocks;
	auto addBlocks = [&](DomainSignature &ds, Side main, Side s) {
		std::bitset<4> neumann = ds.neumannRelative(main);
		int            j       = ds.index(main);
		bool           zp      = ds.zero_patch;
		if (ds.hasFineNbr(main + s)) {
			MatrixBlock b(ds.index(main + s), j, main, s, neumann, zp, BlockType::coarse);
			MatrixBlock c(ds.indexRefinedLeft(main + s), j, main, s, neumann, zp,
			              BlockType::coarse_out_left);
			MatrixBlock d(ds.indexRefinedRight(main + s), j, main, s, neumann, zp,
			              BlockType::coarse_out_right);
			blocks.insert(b);
			blocks.insert(c);
			blocks.insert(d);
		} else if (ds.hasCoarseNbr(main + s)) {
			MatrixBlock b(ds.index(main + s), j, main, s, neumann, zp, BlockType::fine);
			blocks.insert(b);
			if (ds.leftOfCoarse(main + s)) {
				MatrixBlock c(ds.indexCenter(main + s), j, main, s, neumann, zp,
				              BlockType::fine_out_left);
				blocks.insert(c);
			} else {
				MatrixBlock c(ds.indexCenter(main + s), j, main, s, neumann, zp,
				              BlockType::fine_out_right);
				blocks.insert(c);
			}
		} else {
			MatrixBlock b(ds.index(main + s), j, main, s, neumann, zp, BlockType::plain);
			blocks.insert(b);
		}
	};
	// left block
	{
		Side iface_s = Side::north;
		if (y_axis) {
			iface_s = Side::east;
		}
		if (left.index(iface_s) == id_local) {
			// north block
			{
				Side s = Side::north;
				addBlocks(left, iface_s, s);
			}
			// east block
			{
				Side s = Side::east;
				if (left.hasNbr(iface_s + s)) {
					addBlocks(left, iface_s, s);
				}
			}
			// south block
			{
				Side s = Side::south;
				if (left.hasNbr(iface_s + s)) {
					addBlocks(left, iface_s, s);
				}
			}
			// west block
			{
				Side s = Side::west;
				if (left.hasNbr(iface_s + s)) {
					addBlocks(left, iface_s, s);
				}
			}
		}
	}
	// right block
	{
		Side iface_s = Side::south;
		if (y_axis) {
			iface_s = Side::west;
		}
		if (right.index(iface_s) == id_local) {
			// north block
			{
				Side           s = Side::north;
				std::bitset<4> neumann;
				addBlocks(right, iface_s, s);
			}
			// east block
			{
				Side s = Side::east;
				if (right.hasNbr(iface_s + s)) {
					addBlocks(right, iface_s, s);
				}
			}
			// south block
			{
				Side s = Side::south;
				if (right.hasNbr(iface_s + s)) {
					addBlocks(right, iface_s, s);
				}
			}
			// west block
			{
				Side s = Side::west;
				if (right.hasNbr(iface_s + s)) {
					addBlocks(right, iface_s, s);
				}
			}
		}
	}

	return blocks;
}
set<int> Iface::getPins()
{
	set<MatrixBlock> blocks = getGlobalColBlocks();
	set<int>         pins;
	for (MatrixBlock b : blocks) {
		pins.insert(b.i);
	}
	return pins;
}
