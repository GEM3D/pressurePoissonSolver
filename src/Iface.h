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

#ifndef IFACE_H
#define IFACE_H
#include "BufferWriter.h"
#include "IfaceType.h"
#include "Side.h"
#include <bitset>
#include <map>
#include <set>
#include <vector>
template <size_t D> struct Iface {
	IfaceType                           type;
	Side<D>                             s;
	std::array<int, Side<D>::num_sides> ids;
	std::array<int, Side<D>::num_sides> local_id;
	std::array<int, Side<D>::num_sides> global_id;
	std::bitset<Side<D>::num_sides>     neumann;
	Iface()
	{
		ids.fill(-1);
		local_id.fill(-1);
		global_id.fill(-1);
	}
	Iface(std::array<int, Side<D>::num_sides> ids, IfaceType type, Side<D> s,
	      std::bitset<Side<D>::num_sides> neumann)
	{
		this->ids = ids;
		local_id.fill(-1);
		global_id.fill(-1);
		this->type    = type;
		this->s       = s;
		this->neumann = neumann;
	}
};
template <size_t D> struct IfaceSet : public Serializable {
	int                   id        = -1;
	int                   id_local  = -1;
	int                   id_global = -1;
	std::vector<Iface<D>> ifaces;
	std::set<int>         getNbrs() const
	{
		std::set<int> retval;
		for (const Iface<D> &iface : ifaces) {
			for (const int i : iface.ids) {
				if (i != -1 && i != id) { retval.insert(i); }
			}
		}
		return retval;
	}
	void insert(Iface<D> i)
	{
		ifaces.push_back(i);
	}
	void insert(IfaceSet<D> ifs)
	{
		id = ifs.id;
		for (Iface<D> &i : ifs.ifaces) {
			ifaces.push_back(i);
		}
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		id_local = rev_map.at(id);
		for (Iface<D> &iface : ifaces) {
			for (int i = 0; i < Side<D>::num_sides; i++) {
				if (iface.ids[i] != -1) { iface.local_id[i] = rev_map.at(iface.ids[i]); }
			}
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		id_global = rev_map.at(id_local);
		for (Iface<D> &iface : ifaces) {
			for (int i = 0; i < Side<D>::num_sides; i++) {
				if (iface.local_id[i] != -1) { iface.global_id[i] = rev_map.at(iface.local_id[i]); }
			}
		}
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << id;
		writer << id_global;
		int size = ifaces.size();
		writer << size;
		for (const Iface<D> &i : ifaces) {
			writer << i;
		}
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> id;
		reader >> id_global;
		int size = 0;
		reader >> size;
		for (int i = 0; i < size; i++) {
			Iface<D> iface;
			reader >> iface;
			ifaces.push_back(iface);
		}
		return reader.getPos();
	}
};
#endif
