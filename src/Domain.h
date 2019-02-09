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

#ifndef DOMAIN_H
#define DOMAIN_H
#include "BufferWriter.h"
#include "Serializable.h"
#include "Side.h"
#include <array>
#include <bitset>
#include <map>
#include <memory>
#include <vector>
enum class NbrType { Normal, Coarse, Fine };
template <size_t D> struct NbrInfo;
template <size_t D> struct NormalNbrInfo;
template <size_t D> struct CoarseNbrInfo;
template <size_t D> struct FineNbrInfo;
template <size_t D> struct Domain : public Serializable {
	/**
	 * @brief The domain's own id
	 */
	int id        = 0;
	int id_local  = 0;
	int id_global = 0;
	int n         = 10;

	int refine_level  = 1;
	int parent_id     = -1;
	int oct_on_parent = -1;

	std::array<int, Orthant<D>::num_orthants> child_id;

	std::bitset<2 * D> neumann;
	bool               zero_patch = false;

	std::array<double, D> starts;
	std::array<double, D> lengths;

	std::array<NbrInfo<D> *, Side<D>::num_sides> nbr_info;

	Domain()
	{
		starts.fill(0);
		lengths.fill(1);
		nbr_info.fill(nullptr);
		child_id.fill(-1);
	}
	~Domain();
	friend bool operator<(const Domain &l, const Domain &r)
	{
		return l.id < r.id;
	}
	NbrInfo<D> *&     getNbrInfoPtr(Side<D> s);
	NbrType           getNbrType(Side<D> s) const;
	NormalNbrInfo<D> &getNormalNbrInfo(Side<D> s) const;
	CoarseNbrInfo<D> &getCoarseNbrInfo(Side<D> s) const;
	FineNbrInfo<D> &  getFineNbrInfo(Side<D> s) const;
	inline bool       hasNbr(Side<D> s) const;
	inline bool       isNeumann(Side<D> s) const;
	void              setLocalNeighborIndexes(std::map<int, int> &rev_map);
	void              setGlobalNeighborIndexes(std::map<int, int> &rev_map);
	void              setNeumann();
	std::vector<int>  getNbrIds();
	int               serialize(char *buffer) const;
	int               deserialize(char *buffer);
	void              setPtrs(std::map<int, std::shared_ptr<Domain>> &domains);
	void              updateRank(int rank);
	inline bool       hasChildren()
	{
		return child_id[0] != -1;
	}
};
template <size_t D> class NbrInfo : virtual public Serializable
{
	public:
	virtual ~NbrInfo()                                                          = default;
	virtual NbrType getNbrType()                                                = 0;
	virtual void    getNbrIds(std::vector<int> &nbr_ids)                        = 0;
	virtual void    setGlobalIndexes(std::map<int, int> &rev_map)               = 0;
	virtual void    setLocalIndexes(std::map<int, int> &rev_map)                = 0;
	virtual void    setPtrs(std::map<int, std::shared_ptr<Domain<D>>> &domains) = 0;
	virtual void    updateRankOnNeighbors(int new_rank, Side<D> s)              = 0;
};
template <size_t D> class NormalNbrInfo : public NbrInfo<D>
{
	public:
	std::shared_ptr<Domain<D>> ptr          = nullptr;
	int                        rank         = 0;
	int                        id           = 0;
	int                        local_index  = 0;
	int                        global_index = 0;
	NormalNbrInfo() {}
	~NormalNbrInfo() = default;
	NormalNbrInfo(int id)
	{
		this->id = id;
	}
	NbrType getNbrType()
	{
		return NbrType::Normal;
	}
	void getNbrIds(std::vector<int> &nbr_ids)
	{
		nbr_ids.push_back(id);
	};
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		global_index = rev_map.at(local_index);
	}
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		local_index = rev_map.at(id);
	}
	void setPtrs(std::map<int, std::shared_ptr<Domain<D>>> &domains)
	{
		try {
			ptr = domains.at(id);
		} catch (std::out_of_range) {
			ptr = nullptr;
		}
	}
	void updateRank(int new_rank)
	{
		rank = new_rank;
	}
	void updateRankOnNeighbors(int new_rank, Side<D> s)
	{
		ptr->getNormalNbrInfo(s.opposite()).updateRank(new_rank);
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << rank;
		writer << id;
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> rank;
		reader >> id;
		return reader.getPos();
	}
};
template <size_t D> class CoarseNbrInfo : public NbrInfo<D>
{
	public:
	std::shared_ptr<Domain<D>> ptr;
	int                        rank = 0;
	int                        id;
	int                        local_index;
	int                        global_index;
	int                        quad_on_coarse;
	CoarseNbrInfo()  = default;
	~CoarseNbrInfo() = default;
	CoarseNbrInfo(int id, int quad_on_coarse)
	{
		this->id             = id;
		this->quad_on_coarse = quad_on_coarse;
	}
	NbrType getNbrType()
	{
		return NbrType::Coarse;
	}
	void getNbrIds(std::vector<int> &nbr_ids)
	{
		nbr_ids.push_back(id);
	};
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		global_index = rev_map.at(local_index);
	}
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		local_index = rev_map.at(id);
	}
	void setPtrs(std::map<int, std::shared_ptr<Domain<D>>> &domains)
	{
		try {
			ptr = domains.at(id);
		} catch (std::out_of_range) {
			ptr = nullptr;
		}
	}
	void updateRank(int new_rank)
	{
		rank = new_rank;
	}
	void updateRankOnNeighbors(int new_rank, Side<D> s);
	int  serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << rank;
		writer << id;
		writer << quad_on_coarse;
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> rank;
		reader >> id;
		reader >> quad_on_coarse;
		return reader.getPos();
	}
};
template <size_t D> class FineNbrInfo : public NbrInfo<D>
{
	public:
	std::array<std::shared_ptr<Domain<D>>, Orthant<D - 1>::num_orthants> ptrs;
	std::array<int, Orthant<D - 1>::num_orthants>                        ranks;
	std::array<int, Orthant<D - 1>::num_orthants>                        ids;
	std::array<int, Orthant<D - 1>::num_orthants>                        global_indexes;
	std::array<int, Orthant<D - 1>::num_orthants>                        local_indexes;
	FineNbrInfo()
	{
		ptrs.fill(nullptr);
		ranks.fill(0);
	}
	~FineNbrInfo() = default;
	FineNbrInfo(std::array<int, Orthant<D>::num_orthants / 2> ids)
	{
		this->ids = ids;
	}
	NbrType getNbrType()
	{
		return NbrType::Fine;
	}
	void getNbrIds(std::vector<int> &nbr_ids)
	{
		for (size_t i = 0; i < ids.size(); i++) {
			nbr_ids.push_back(ids[i]);
		}
	};
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		for (size_t i = 0; i < global_indexes.size(); i++) {
			global_indexes[i] = rev_map.at(local_indexes[i]);
		}
	}
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		for (size_t i = 0; i < local_indexes.size(); i++) {
			local_indexes[i] = rev_map.at(ids[i]);
		}
	}
	void setPtrs(std::map<int, std::shared_ptr<Domain<D>>> &domains)
	{
		for (size_t i = 0; i < ids.size(); i++) {
			try {
				ptrs[i] = domains.at(ids[i]);
			} catch (std::out_of_range) {
				ptrs[i] = nullptr;
			}
		}
	}
	void updateRank(int new_rank, int quad_on_coarse)
	{
		ranks[quad_on_coarse] = new_rank;
	}

	void updateRankOnNeighbors(int new_rank, Side<D> s)
	{
		for (size_t i = 0; i < ptrs.size(); i++) {
			ptrs[i]->getCoarseNbrInfo(s.opposite()).updateRank(new_rank);
		}
	}
	int serialize(char *buffer) const
	{
		BufferWriter writer(buffer);
		writer << ranks;
		writer << ids;
		return writer.getPos();
	}
	int deserialize(char *buffer)
	{
		BufferReader reader(buffer);
		reader >> ranks;
		reader >> ids;
		return reader.getPos();
	}
};
template <size_t D> inline void CoarseNbrInfo<D>::updateRankOnNeighbors(int new_rank, Side<D> s)
{
	ptr->getFineNbrInfo(s.opposite()).updateRank(new_rank, quad_on_coarse);
}
template <size_t D> inline Domain<D>::~Domain()
{
	/*
	for (NbrInfo *info : nbr_info) {
	    delete info;
	}
	*/
}
template <size_t D> inline NbrInfo<D> *&Domain<D>::getNbrInfoPtr(Side<D> s)
{
	return nbr_info[s.toInt()];
}
template <size_t D> inline NbrType Domain<D>::getNbrType(Side<D> s) const
{
	return nbr_info[s.toInt()]->getNbrType();
}
template <size_t D> inline NormalNbrInfo<D> &Domain<D>::getNormalNbrInfo(Side<D> s) const
{
	return *(NormalNbrInfo<D> *) nbr_info[s.toInt()];
}
template <size_t D> inline CoarseNbrInfo<D> &Domain<D>::getCoarseNbrInfo(Side<D> s) const
{
	return *(CoarseNbrInfo<D> *) nbr_info[s.toInt()];
}
template <size_t D> inline FineNbrInfo<D> &Domain<D>::getFineNbrInfo(Side<D> s) const
{
	return *(FineNbrInfo<D> *) nbr_info[s.toInt()];
}
template <size_t D> inline bool Domain<D>::hasNbr(Side<D> s) const
{
	return nbr_info[s.toInt()] != nullptr;
}
template <size_t D> inline bool Domain<D>::isNeumann(Side<D> s) const
{
	return neumann[s.toInt()];
}
template <size_t D> inline void Domain<D>::setLocalNeighborIndexes(std::map<int, int> &rev_map)
{
	id_local = rev_map.at(id);
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setLocalIndexes(rev_map); }
	}
}
template <size_t D> inline void Domain<D>::setGlobalNeighborIndexes(std::map<int, int> &rev_map)
{
	id_global = rev_map.at(id_local);
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setGlobalIndexes(rev_map); }
	}
}
template <size_t D> inline void Domain<D>::setNeumann()
{
	for (size_t q = 0; q < neumann.size(); q++) {
		neumann[q] = !hasNbr(static_cast<Side<D>>(q));
	}
}
template <size_t D> inline std::vector<int> Domain<D>::getNbrIds()
{
	std::vector<int> retval;
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->getNbrIds(retval); }
	}
	return retval;
}

template <size_t D> inline int Domain<D>::serialize(char *buffer) const
{
	BufferWriter writer(buffer);
	writer << id;
	writer << n;
	writer << refine_level;
	writer << parent_id;
	writer << oct_on_parent;
	writer << child_id;
	writer << neumann;
	writer << zero_patch;
	writer << starts;
	writer << lengths;
	std::bitset<Side<D>::num_sides> has_nbr;
	for (size_t i = 0; i < Side<D>::num_sides; i++) {
		has_nbr[i] = nbr_info[i] != nullptr;
	}
	writer << has_nbr;
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) {
			NbrType type = getNbrType(s);
			writer << type;
			switch (type) {
				case NbrType::Normal: {
					NormalNbrInfo<D> tmp = getNormalNbrInfo(s);
					writer << tmp;
				} break;
				case NbrType::Fine: {
					FineNbrInfo<D> tmp = getFineNbrInfo(s);
					writer << tmp;
				} break;
				case NbrType::Coarse: {
					CoarseNbrInfo<D> tmp = getCoarseNbrInfo(s);
					writer << tmp;
				} break;
			}
		}
	}
	return writer.getPos();
}
template <size_t D> inline int Domain<D>::deserialize(char *buffer)
{
	BufferReader reader(buffer);
	reader >> id;
	reader >> n;
	reader >> refine_level;
	reader >> parent_id;
	reader >> oct_on_parent;
	reader >> child_id;
	reader >> neumann;
	reader >> zero_patch;
	reader >> starts;
	reader >> lengths;
	std::bitset<Side<D>::num_sides> has_nbr;
	reader >> has_nbr;
	for (size_t i = 0; i < Side<D>::num_sides; i++) {
		if (has_nbr[i]) {
			NbrType type;
			reader >> type;
			NbrInfo<D> *info = nullptr;
			switch (type) {
				case NbrType::Normal:
					info = new NormalNbrInfo<D>();
					reader >> *(NormalNbrInfo<D> *) info;
					break;
				case NbrType::Fine:
					info = new FineNbrInfo<D>();
					reader >> *(FineNbrInfo<D> *) info;
					break;
				case NbrType::Coarse:
					info = new CoarseNbrInfo<D>();
					reader >> *(CoarseNbrInfo<D> *) info;
					break;
			}
			nbr_info[i] = info;
		}
	}
	return reader.getPos();
}
template <size_t D> inline void Domain<D>::setPtrs(std::map<int, std::shared_ptr<Domain>> &domains)
{
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setPtrs(domains); }
	}
}
template <size_t D> inline void Domain<D>::updateRank(int rank)
{
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->updateRankOnNeighbors(rank, s); }
	}
}
#endif
