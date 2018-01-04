#ifndef IFACE_H
#define IFACE_H
#include "IfaceType.h"
#include "Side.h"
#include <map>
#include <set>
#include <vector>
struct Iface {
	IfaceType type;
	Side      s;
	std::array<int, 6> g_id      = {{-1, -1, -1, -1, -1, -1}};
	std::array<int, 6> local_id  = {{-1, -1, -1, -1, -1, -1}};
	std::array<int, 6> global_id = {{-1, -1, -1, -1, -1, -1}};
};
struct IfaceSet {
	int                gid       = -1;
	int                id_local  = -1;
	int                id_global = -1;
	std::vector<Iface> ifaces;
	std::set<int>      getNbrs() const
	{
		std::set<int> retval;
		for (const Iface &iface : ifaces) {
			for (const int i : iface.g_id) {
				if (i != -1 && i != gid) {
					retval.insert(i);
				}
			}
		}
		return retval;
	}
	void insert(Iface i) { ifaces.push_back(i); }
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		id_local = rev_map.at(gid);
		for (Iface &iface : ifaces) {
			for (int i = 0; i < 6; i++) {
				if (iface.g_id[i] != -1) {
					iface.local_id[i] = rev_map.at(iface.g_id[i]);
				}
			}
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		id_global = rev_map.at(id_local);
		for (Iface &iface : ifaces) {
			for (int i = 0; i < 6; i++) {
				if (iface.local_id[i] != -1) {
					iface.global_id[i] = rev_map.at(iface.local_id[i]);
				}
			}
		}
	}
};
#endif
