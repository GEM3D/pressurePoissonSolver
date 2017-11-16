#include "DomainSignatureCollection.h"
#include "MyTypeDefs.h"
#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <utility>
using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
DomainSignatureCollection::DomainSignatureCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm,
                                                     string file_name, int rank)
{
	this->comm = comm;
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
DomainSignatureCollection::DomainSignatureCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm,
                                                     int d_x, int d_y, int rank)
{
	this->comm         = comm;
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
			int index = ds.gid(s);
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
	indexIfacesGlobal();
	indexDomainIfacesLocal();
}
DomainSignatureCollection::DomainSignatureCollection(Teuchos::RCP<const Teuchos::Comm<int>> comm,
                                                     int d_x, int d_y, int rank, bool amr)
{
	this->comm         = comm;
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
				if (d.hasNbr(s) && d.gid(s) == -1) {
					// a new edge that we have not assigned an index to
					d.gid(s) = curr_i;
					curr_i++;

					// fine case
					if (d.hasFineNbr(s)) {
						DomainSignature &nbr_left  = domains.at(d.nbr(s));
						DomainSignature &nbr_right = domains.at(d.nbrRight(s));

						// set center indexes
						nbr_left.gidCenter(!s)  = d.gid(s);
						nbr_right.gidCenter(!s) = d.gid(s);

						// set left and right indexes index
						nbr_left.gid(!s) = curr_i;
						curr_i++;
						nbr_right.gid(!s) = curr_i;
						curr_i++;

						// set refined indexes
						d.gidRefinedLeft(s)  = nbr_left.gid(!s);
						d.gidRefinedRight(s) = nbr_right.gid(!s);

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

							nbr.gidRefinedLeft(!s) = d.gid(s);

							nbr.gidRefinedRight(!s) = curr_i;
							buddy.gid(s)            = curr_i;
							curr_i++;

							d.gidCenter(s)     = curr_i;
							nbr.gid(!s)        = curr_i;
							buddy.gidCenter(s) = curr_i;
							curr_i++;
						} else {
							DomainSignature &buddy = domains.at(nbr.nbr(!s));
							buddy_id               = buddy.id;

							nbr.gidRefinedRight(!s) = d.gid(s);

							nbr.gidRefinedLeft(!s) = curr_i;
							buddy.gid(s)           = curr_i;
							curr_i++;

							d.gidCenter(s)     = curr_i;
							nbr.gid(!s)        = curr_i;
							buddy.gidCenter(s) = curr_i;
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
						nbr.gid(!s)          = d.gid(s);
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
	zoltanBalanceDomains();
	indexDomainIfacesLocal();
	zoltanBalanceIfaces();
	indexIfacesLocal();
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
	ostringstream oss;
	oss << "I have " << ifaces.size() << " ifaces: ";

#if DD_DEBUG
	int  prev  = -100;
	bool range = false;
	for (auto &p : ifaces) {
		int curr = p.second.id;
		if (curr != prev + 1 && !range) {
			oss << curr << "-";
			range = true;
		} else if (curr != prev + 1 && range) {
			oss << prev << " " << curr << "-";
		}
		prev = curr;
	}

	oss << prev << "\n";
#endif
	oss << endl;
	cerr << oss.str();
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

void DomainSignatureCollection::indexIfacesGlobal()
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
			Iface &iface        = ifaces[i];
			rev_map[i]          = curr_i;
			ifaces[i].id_global = curr_i;
			curr_i++;
			set<MatrixBlock> blocks = iface.getGidRowBlocks();
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
			p.second.setGlobalIndexes(rev_map);
		}
	}
	iface_rev_map          = rev_map;
	iface_map_vec          = map_vec;
	iface_off_proc_map_vec = off_proc_map_vec;
}
void DomainSignatureCollection::indexIfacesLocal()
{
	int         curr_i = 0;
	vector<int> map_vec;
	vector<int> off_proc_map_vec;
	vector<int> off_proc_map_vec_send;
	map<int, int> rev_map;
	if (!ifaces.empty()) {
		set<int> todo;
		for (auto &p : ifaces) {
			todo.insert(p.first);
		}
		while (!todo.empty()) {
			deque<int> queue;
			set<int>   enqueued;
			queue.push_back(*todo.begin());
			enqueued.insert(*todo.begin());
			deque<int> off_proc_ifaces;
			while (!queue.empty()) {
				int i = queue.front();
				todo.erase(i);
				queue.pop_front();
				map_vec.push_back(i);
				Iface &iface       = ifaces[i];
				rev_map[i]         = curr_i;
				ifaces[i].id_local = curr_i;
				curr_i++;
				set<MatrixBlock> blocks = iface.getGidRowBlocks();
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
		}
	}
	// sort off proc indeces by proc that they reside on
	// creat tpetra map get procs that indeces reside on
	map_type iface_map = map_type(-1, &map_vec[0], map_vec.size(), 0, comm);

	set<int> neighbors;
	map<int, set<int>> proc_recv;
	// RECIEVING
	{
		vector<int>             procs(off_proc_map_vec.size());
		Teuchos::ArrayView<int> inds_view(&off_proc_map_vec[0], off_proc_map_vec.size());
		Teuchos::ArrayView<int> procs_view(&procs[0], procs.size());
		iface_map.getRemoteIndexList(inds_view, procs_view);

		// sort
		vector<pair<int, int>> proc_ind(procs.size());
		for (size_t i = 0; i < off_proc_map_vec.size(); i++) {
			proc_ind[i] = make_pair(procs[i], off_proc_map_vec[i]);
		}

		sort(proc_ind.begin(), proc_ind.end());

		for (size_t i = 0; i < proc_ind.size(); i++) {
			procs[i]            = proc_ind[i].first;
			off_proc_map_vec[i] = proc_ind[i].second;
			neighbors.insert(procs[i]);
			proc_recv[procs[i]].insert(off_proc_map_vec[i]);
		}
		iface_off_proc_vec = procs;
	}

	// map off proc
	for (int i : off_proc_map_vec) {
		// TODO map in increasing order of neighboring proc
		rev_map[i] = curr_i;
		curr_i++;
	}
	for (auto &p : ifaces) {
		p.second.setLocalIndexes(rev_map);
	}
	// SENDING
	vector<int> src;
	vector<int> dest;
	set<int>    uniq_dest;
	for (auto &p : ifaces) {
		set<MatrixBlock> blocks = p.second.getGlobalColBlocks();
		for (MatrixBlock b : blocks) {
			if (ifaces.count(b.i) == 0) {
				src.push_back(b.j);
				dest.push_back(b.i);
				uniq_dest.insert(b.i);
			}
		}
	}
	// get procs
	map<int, set<int>> proc_send;
	{
		vector<int> inds;
		for (int i : uniq_dest) {
			inds.push_back(i);
		}
		vector<int>             procs(inds.size());
		Teuchos::ArrayView<int> inds_view(&inds[0], inds.size());
		Teuchos::ArrayView<int> procs_view(&procs[0], procs.size());
		iface_map.getRemoteIndexList(inds_view, procs_view);

		map<int, int> ind_proc;
		for (size_t i = 0; i < procs.size(); i++) {
			ind_proc[inds[i]] = procs[i];
		}
		// separate into procs
		for (size_t i = 0; i < dest.size(); i++) {
			int d = dest[i];
			int p = ind_proc[d];
			proc_send[p].insert(src[i]);
			neighbors.insert(p);
		}
	}

	// create amgxmap
	amgxmap = AmgxMap(neighbors.size());
	map<int, int> neighbors_rev;
	int n_i = 0;
	for (int proc : neighbors) {
		amgxmap.neighbors[n_i] = proc;
		neighbors_rev[proc]    = n_i;
		n_i++;
	}
	// create send maps
	for (auto &p : proc_send) {
		int n_i                 = neighbors_rev[p.first];
		int size                = p.second.size();
		amgxmap.send_sizes[n_i] = size;
		ArrayRCP<int> lmap(size);
		int           i = 0;
		for (int ind : p.second) {
			lmap[i] = rev_map.at(ind);
			i++;
		}
		amgxmap.send_maps[n_i] = &lmap[0];
		amgxmap.arrays.insert(lmap);
	}
	// create recv maps
	for (auto &p : proc_recv) {
		int n_i                 = neighbors_rev[p.first];
		int size                = p.second.size();
		amgxmap.recv_sizes[n_i] = size;
		ArrayRCP<int> lmap(size);
		int           i = 0;
		for (int ind : p.second) {
			lmap[i] = rev_map.at(ind);
			i++;
		}
		amgxmap.recv_maps[n_i] = &lmap[0];
		amgxmap.arrays.insert(lmap);
	}
	iface_rev_map               = rev_map;
	iface_map_vec               = map_vec;
	iface_off_proc_map_vec      = off_proc_map_vec;
	iface_off_proc_map_vec_send = off_proc_map_vec_send;
}
AmgxMap::AmgxMap(const AmgxMap &orig, int n) : AmgxMap(orig.num_neighbors)
{
	neighbors = orig.neighbors;
	for (size_t i = 0; i < orig.recv_sizes.size(); i++) {
		recv_sizes[i] = orig.recv_sizes[i] * n;
	}
	for (size_t i = 0; i < orig.send_sizes.size(); i++) {
		send_sizes[i] = orig.send_sizes[i] * n;
	}
	for (int i = 0; i < num_neighbors; i++) {
		int size = orig.recv_sizes[i];
		if (size != 0) {
			ArrayRCP<int> lmap(size * n);
			for (int j = 0; j < size; j++) {
				int ind = orig.recv_maps[i][j];
				for (int k = 0; k < n; k++) {
					lmap[j * n + k] = ind * n + k;
				}
			}
			recv_maps[i] = &lmap[0];
			arrays.insert(lmap);
		}
        cerr << endl;
		size = orig.send_sizes[i];
		if (size != 0) {
			ArrayRCP<int> lmap(size * n);
			for (int j = 0; j < size; j++) {
				int ind = orig.send_maps[i][j];
				for (int k = 0; k < n; k++) {
					lmap[j * n + k] = ind * n + k;
				}
			}
			send_maps[i] = &lmap[0];
			arrays.insert(lmap);
		}
        cerr << endl;
	}
}
void DomainSignatureCollection::indexDomainIfacesLocal()
{
	vector<int> map_vec;
	map<int, int> rev_map;
	if (!domains.empty()) {
		int      curr_i = 0;
		set<int> todo;
		for (auto &p : domains) {
			todo.insert(p.first);
		}
		while (!todo.empty()) {
			deque<int> queue;
			set<int>   enqueued;
			queue.push_back(*todo.begin());
			enqueued.insert(*todo.begin());
			deque<int> off_proc_ifaces;
			while (!queue.empty()) {
				int i = queue.front();
				queue.pop_front();
				todo.erase(i);
				DomainSignature &ds = domains[i];
				for (int g_id : ds.g_id) {
					if (g_id != -1 && rev_map.count(g_id) == 0) {
						rev_map[g_id] = curr_i;
						map_vec.push_back(g_id);
						curr_i++;
					}
				}
				for (int g_id : ds.g_id_refined) {
					if (g_id != -1 && rev_map.count(g_id) == 0) {
						rev_map[g_id] = curr_i;
						map_vec.push_back(g_id);
						curr_i++;
					}
				}
				for (int g_id : ds.g_id_center) {
					if (g_id != -1 && rev_map.count(g_id) == 0) {
						rev_map[g_id] = curr_i;
						map_vec.push_back(g_id);
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
		}
		for (auto &p : domains) {
			p.second.setLocalIndexes(rev_map);
		}
	}
	domain_rev_map = rev_map;
	domain_map_vec = map_vec;
}
std::set<MatrixBlock> Iface::getRowBlocks(int *id, DsMemPtr normal)
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
	int i = *id;
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
			int            j = (left.*normal)(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, false, getBlockType(true));
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = (left.*normal)(main);
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
			int  j    = (left.*normal)(main);
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
			int  j    = (left.*normal)(main);
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
			int            j = (right.*normal)(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, false, getBlockType(false));
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = (right.*normal)(main);
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
			int  j    = (right.*normal)(main);
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
			int  j    = (right.*normal)(main);
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
			int            j = (extra.*normal)(main);
			MatrixBlock    b(i, j, iface_s, Side::north, neumann, false, BlockType::fine_out_right);
			blocks.insert(b);
		}
		// east block
		{
			Side s    = Side::west;
			Side rel  = Side::east;
			Side main = iface_s + rel;
			int  j    = (extra.*normal)(main);
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
			int  j    = (extra.*normal)(main);
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
			int  j    = (extra.*normal)(main);
			if (j != -1) {
				std::bitset<4> neumann;
				MatrixBlock    b(i, j, main, s, neumann, false, BlockType::fine_out_right);
				blocks.insert(b);
			}
		}
	}
	return blocks;
}
std::set<MatrixBlock> Iface::getGidRowBlocks() { return getRowBlocks(&id, &DomainSignature::gid); }
std::set<MatrixBlock> Iface::getGlobalRowBlocks()
{
	return getRowBlocks(&id_global, &DomainSignature::globalIndex);
}
std::set<MatrixBlock> Iface::getRowBlocks()
{
	return getRowBlocks(&id_local, &DomainSignature::index);
}
std::set<MatrixBlock> Iface::getGlobalColBlocks()
{
	std::set<MatrixBlock> blocks;
	auto addBlocks = [&](DomainSignature &ds, Side main, Side s) {
		std::bitset<4> neumann = ds.neumannRelative(main);
		int            j       = ds.gid(main);
		bool           zp      = ds.zero_patch;
		if (ds.hasFineNbr(main + s)) {
			MatrixBlock b(ds.gid(main + s), j, main, s, neumann, zp, BlockType::coarse);
			MatrixBlock c(ds.gidRefinedLeft(main + s), j, main, s, neumann, zp,
			              BlockType::coarse_out_left);
			MatrixBlock d(ds.gidRefinedRight(main + s), j, main, s, neumann, zp,
			              BlockType::coarse_out_right);
			blocks.insert(b);
			blocks.insert(c);
			blocks.insert(d);
		} else if (ds.hasCoarseNbr(main + s)) {
			MatrixBlock b(ds.gid(main + s), j, main, s, neumann, zp, BlockType::fine);
			blocks.insert(b);
			if (ds.leftOfCoarse(main + s)) {
				MatrixBlock c(ds.gidCenter(main + s), j, main, s, neumann, zp,
				              BlockType::fine_out_left);
				blocks.insert(c);
			} else {
				MatrixBlock c(ds.gidCenter(main + s), j, main, s, neumann, zp,
				              BlockType::fine_out_right);
				blocks.insert(c);
			}
		} else {
			MatrixBlock b(ds.gid(main + s), j, main, s, neumann, zp, BlockType::plain);
			blocks.insert(b);
		}
	};
	// left block
	{
		Side iface_s = Side::north;
		if (y_axis) {
			iface_s = Side::east;
		}
		if (left.gid(iface_s) == id) {
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
		if (right.gid(iface_s) == id) {
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
set<int> Iface::getPins()
{
	set<MatrixBlock> blocks = getGidRowBlocks();
	set<int>         pins;
	for (MatrixBlock b : blocks) {
		pins.insert(b.j);
	}
	return pins;
}
