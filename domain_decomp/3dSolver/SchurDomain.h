#ifndef SCHURDOMAIN_H
#define SCHURDOMAIN_H
#include "Domain.h"
#include "Iface.h"
#include <deque>
class IfaceInfo
{
	public:
	int          id;
	int          local_index;
	int          global_index;
	bool         own;
	virtual void getIds(std::vector<int> &ids) = 0;
	virtual void getIdsAndTypes(std::deque<int> &ids, std::deque<IfaceType> &types,
	                            std::deque<Side> &sides, std::deque<int> &ranks, Side s)
	= 0;
	virtual void getIdxAndTypes(std::deque<int> &idx, std::deque<IfaceType> &types) = 0;
	virtual void setLocalIndexes(const std::map<int, int> &rev_map)                 = 0;
	virtual void setGlobalIndexes(const std::map<int, int> &rev_map)                = 0;
};
class NormalIfaceInfo : public IfaceInfo
{
	public:
	NormalNbrInfo nbr_info;
	NormalIfaceInfo()
	{
		id           = 0;
		local_index  = 0;
		global_index = 0;
	}
	NormalIfaceInfo(Domain &d, Side s)
	{
		nbr_info = d.getNormalNbrInfo(s);
		own      = d.getNormalNbrInfo(s).ptr != nullptr;
		switch (s.toInt()) {
			case Side::west:
			case Side::south:
			case Side::bottom:
				id = d.id * 6 + s.toInt();
				break;
			case Side::east:
			case Side::north:
			case Side::top:
				id  = d.getNormalNbrInfo(s).id * 6 + s.opposite().toInt();
				own = true;
		}
	}
	void getIds(std::vector<int> &ids)
	{
		ids.push_back(id);
	}
	void getIdsAndTypes(std::deque<int> &ids, std::deque<IfaceType> &types, std::deque<Side> &sides,
	                    std::deque<int> &ranks, Side s)
	{
		ids.push_back(id);
		types.push_back(IfaceType::normal);
		sides.push_back(s);
		ranks.push_back(nbr_info.rank);
	}
	void getIdxAndTypes(std::deque<int> &idx, std::deque<IfaceType> &types)
	{
		idx.push_back(local_index);
		types.push_back(IfaceType::normal);
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		local_index = rev_map.at(id);
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		global_index = rev_map.at(local_index);
	}
};
class CoarseIfaceInfo : public IfaceInfo
{
	public:
	CoarseNbrInfo nbr_info;
	int           quad_on_coarse;
	int           coarse_id;
	int           coarse_local_index;
	int           coarse_global_index;
	CoarseIfaceInfo(Domain &d, Side s)
	{
		nbr_info                = d.getCoarseNbrInfo(s);
		id                      = d.id * 6 + s.toInt();
		CoarseNbrInfo &nbr_info = d.getCoarseNbrInfo(s);
		quad_on_coarse          = nbr_info.quad_on_coarse;
		coarse_id               = nbr_info.id * 6 + s.opposite().toInt();
		own                     = true;
	}
	void getIdsAndTypes(std::deque<int> &ids, std::deque<IfaceType> &types, std::deque<Side> &sides,
	                    std::deque<int> &ranks, Side s)
	{
		ids.push_back(id);
		ids.push_back(coarse_id);
		types.push_back(IfaceType::fine_to_fine_0 + quad_on_coarse);
		types.push_back(IfaceType::fine_to_coarse_0 + quad_on_coarse);
		sides.push_back(s);
		sides.push_back(s);
		ranks.push_back(nbr_info.rank);
		ranks.push_back(nbr_info.rank);
	}
	void getIdxAndTypes(std::deque<int> &idx, std::deque<IfaceType> &types)
	{
		idx.push_back(local_index);
		idx.push_back(coarse_local_index);
		types.push_back(IfaceType::fine_to_fine_0 + quad_on_coarse);
		types.push_back(IfaceType::fine_to_coarse_0 + quad_on_coarse);
	}
	void getIds(std::vector<int> &ids)
	{
		ids.push_back(id);
		ids.push_back(coarse_id);
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		local_index        = rev_map.at(id);
		coarse_local_index = rev_map.at(coarse_id);
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		global_index        = rev_map.at(local_index);
		coarse_global_index = rev_map.at(coarse_local_index);
	}
};
class FineIfaceInfo : public IfaceInfo
{
	public:
	FineNbrInfo        nbr_info;
	std::array<int, 4> fine_ids;
	std::array<int, 4> fine_local_indexes;
	std::array<int, 4> fine_global_indexes;
	FineIfaceInfo(Domain &d, Side s)
	{
		nbr_info              = d.getFineNbrInfo(s);
		id                    = d.id * 6 + s.toInt();
		FineNbrInfo &nbr_info = d.getFineNbrInfo(s);
		for (int i = 0; i < 4; i++) {
			fine_ids[i] = nbr_info.ids[i] * 6 + s.opposite().toInt();
		}
		own = true;
	}
	void getIdsAndTypes(std::deque<int> &ids, std::deque<IfaceType> &types, std::deque<Side> &sides,
	                    std::deque<int> &ranks, Side s)
	{
		ids.push_back(id);
		types.push_back(IfaceType::coarse_to_coarse);
		sides.push_back(s);
		for (int i = 0; i < 4; i++) {
			ids.push_back(fine_ids[i]);
			types.push_back(IfaceType::coarse_to_fine_0 + i);
			sides.push_back(s);
			ranks.push_back(nbr_info.ranks[i]);
		}
	}
	void getIdxAndTypes(std::deque<int> &idx, std::deque<IfaceType> &types)
	{
		idx.push_back(local_index);
		types.push_back(IfaceType::coarse_to_coarse);
		for (int i = 0; i < 4; i++) {
			idx.push_back(fine_local_indexes[i]);
			types.push_back(IfaceType::coarse_to_fine_0 + i);
		}
	}
	void getIds(std::vector<int> &ids)
	{
		ids.push_back(id);
		for (int i = 0; i < 4; i++) {
			ids.push_back(fine_ids[i]);
		}
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		local_index = rev_map.at(id);
		for (int i = 0; i < 4; i++) {
			fine_local_indexes[i] = rev_map.at(fine_ids[i]);
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		global_index = rev_map.at(local_index);
		for (int i = 0; i < 4; i++) {
			fine_global_indexes[i] = rev_map.at(fine_local_indexes[i]);
		}
	}
};
struct SchurDomain {
	Domain                     domain;
	int                        local_index = 0;
	int                        n;
	double                     x_length;
	double                     y_length;
	double                     z_length;
	std::bitset<6>             neumann;
	std::array<IfaceInfo *, 6> iface_info
	= {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}};
	std::vector<int> nbr_ids;
	SchurDomain() = default;
	SchurDomain(Domain &d)
	{
		domain      = d;
		local_index = d.id_local;
		n           = d.n;
		x_length    = d.x_length;
		y_length    = d.y_length;
		z_length    = d.z_length;
		neumann     = d.neumann;

		// create iface objects
		for (Side s : Side::getValues()) {
			if (d.hasNbr(s)) {
				switch (d.getNbrType(s)) {
					case NbrType::Normal:
						getIfaceInfoPtr(s) = new NormalIfaceInfo(d, s);
						break;
					case NbrType::Fine:
						getIfaceInfoPtr(s) = new FineIfaceInfo(d, s);
						break;
					case NbrType::Coarse:
						getIfaceInfoPtr(s) = new CoarseIfaceInfo(d, s);
						break;
				}
			}
		}
	}
	IfaceInfo *&getIfaceInfoPtr(Side s)
	{
		return iface_info[s.toInt()];
	}
	NormalIfaceInfo &getNormalIfaceInfo(Side s)
	{
		return *(NormalIfaceInfo *) iface_info[s.toInt()];
	}
	bool hasNbr(Side s)
	{
		return iface_info[s.toInt()] != nullptr;
	}
	void enumerateIfaces(std::map<int, IfaceSet> &                ifaces,
	                     std::map<int, std::pair<int, IfaceSet>> &off_proc_ifaces)
	{
		std::array<int, 6> ids;
		for (Side s : Side::getValues()) {
			if (hasNbr(s)) {
				ids[s.toInt()] = getIfaceInfoPtr(s)->id;
			} else {
				ids[s.toInt()] = -1;
			}
		}
		std::deque<int>       iface_ids;
		std::deque<IfaceType> iface_types;
		std::deque<Side>      iface_sides;
		std::deque<bool>      iface_own;
		std::deque<int>       iface_ranks;
		for (Side s : Side::getValues()) {
			if (hasNbr(s)) {
				int num_added = iface_ids.size();
				getIfaceInfoPtr(s)->getIdsAndTypes(iface_ids, iface_types, iface_sides, iface_ranks,
				                                   s);
				num_added = iface_ids.size() - num_added;
				for (int i = 0; i < num_added; i++) {
					iface_own.push_back(getIfaceInfoPtr(s)->own);
				}
			}
		}
		for (size_t i = 0; i < iface_ids.size(); i++) {
			int       id   = iface_ids[i];
			IfaceType type = iface_types[i];
			Side      s    = iface_sides[i];
			int       rank = iface_ranks[i];
			if (iface_own[i]) {
				IfaceSet &ifs = ifaces[id];
				ifs.id        = id;
				ifs.insert(Iface(ids, type, s, neumann));
			} else {
				IfaceSet &ifs             = off_proc_ifaces[id].second;
				off_proc_ifaces[id].first = rank;
				ifs.id                    = id;
				ifs.insert(Iface(ids, type, s, neumann));
			}
		}
	}
	std::vector<int> getIds()
	{
		std::vector<int> retval;
		for (Side s : Side::getValues()) {
			if (hasNbr(s)) { getIfaceInfoPtr(s)->getIds(retval); }
		}
		return retval;
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		for (Side s : Side::getValues()) {
			if (hasNbr(s)) { getIfaceInfoPtr(s)->setLocalIndexes(rev_map); }
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		for (Side s : Side::getValues()) {
			if (hasNbr(s)) { getIfaceInfoPtr(s)->setGlobalIndexes(rev_map); }
		}
	}
	bool isNeumann(Side s)
	{
		return neumann[s.toInt()];
	}
	int getIfaceLocalIndex(Side s)
	{
		return iface_info[s.toInt()]->local_index;
	}
};
#endif
