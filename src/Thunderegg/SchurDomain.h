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

#ifndef SCHURDOMAIN_H
#define SCHURDOMAIN_H
#include "Domain.h"
#include "Iface.h"
#include <deque>
/**
 * @brief The IfaceInfo class represents the information for an interface on a given side of the
 * patch.
 *
 * The information contained  will be the globally unique ID and the local and global index(es) in
 * the interface vector.
 *
 * @tparam D the number of Cartesian dimensions in the patches.
 */
template <size_t D> class IfaceInfo
{
	public:
	/**
	 * @brief The globally unique ID of the interface.
	 */
	int id;
	/**
	 * @brief the local index in the interface vector.
	 */
	int local_index;
	/**
	 * @brief the global index in the interface vector.
	 */
	int global_index;
	/**
	 * @brief Whether the global location of this particular interface resides on the same
	 * processor.
	 */
	bool own;
	/**
	 * @brief add to a deque of globally unique ids
	 *
	 * @param ids adds to the deque of ids
	 */
	virtual void getIds(std::deque<int> &ids) = 0;
	/**
	 * @brief add to a deque of local interface indexes
	 *
	 * @param ids adds to the deque of local interface indexes
	 */
	virtual void getLocalIndexes(std::deque<int> &idx) = 0;
	/**
	 * @brief add to a deque of global interface indexes
	 *
	 * @param ids adds to the deque of global interface indexes
	 */
	virtual void getGlobalIndexes(std::deque<int> &idx) = 0;
	/**
	 * @brief add to a deque of IfaceTypes
	 *
	 * @param ids adds to the deque of IfaceTypes
	 */
	virtual void getIfaceTypes(std::deque<IfaceType> &types) = 0;
	/**
	 * @brief add to a deque of interface ranks
	 *
	 * @param ids adds to the deque of interface ranks
	 */
	virtual void getRanks(std::deque<int> &ranks) = 0;
	/**
	 * @brief add to a deque of interface ownership
	 *
	 * @param ids adds to the deque of interface ownership
	 */
	virtual void getOwnership(std::deque<bool> &own)             = 0;
	virtual void getIncomingProcs(std::set<int> &incoming_procs) = 0;
	virtual void setLocalIndexes(const std::map<int, int> &rev_map)  = 0;
	virtual void setGlobalIndexes(const std::map<int, int> &rev_map) = 0;
};
template <size_t D> class NormalIfaceInfo : public IfaceInfo<D>
{
	public:
	NormalNbrInfo<D> nbr_info;
	NormalIfaceInfo()
	{
		this->id           = 0;
		this->local_index  = 0;
		this->global_index = 0;
	}
	NormalIfaceInfo(Domain<D> &d, Side<D> s)
	{
		nbr_info = d.getNormalNbrInfo(s);
		if (s.isLowerOnAxis()) {
			this->id = d.id * Side<D>::num_sides + s.toInt();
		} else {
			this->id = d.getNormalNbrInfo(s).id * Side<D>::num_sides + s.opposite().toInt();
		}
		this->own = ((s.toInt() & 0x1) || nbr_info.ptr != nullptr);
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
	}
	void getIfaceTypes(std::deque<IfaceType> &types)
	{
		types.push_back(IfaceType::normal);
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(nbr_info.rank);
	}
	void getOwnership(std::deque<bool> &own)
	{
		own.push_back(this->own);
	}
	void getIncomingProcs(std::set<int> &incoming_procs)
	{
		if (this->own && nbr_info.ptr == nullptr) { incoming_procs.insert(nbr_info.rank); }
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index = rev_map.at(this->id);
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index = rev_map.at(this->local_index);
	}
};
template <size_t D> class CoarseIfaceInfo : public IfaceInfo<D>
{
	public:
	CoarseNbrInfo<D> nbr_info;
	int              quad_on_coarse;
	int              coarse_id;
	int              coarse_local_index;
	int              coarse_global_index;
	bool             coarse_own;
	CoarseIfaceInfo(Domain<D> &d, Side<D> s)
	{
		nbr_info       = d.getCoarseNbrInfo(s);
		this->id       = d.id * Side<D>::num_sides + s.toInt();
		nbr_info       = d.getCoarseNbrInfo(s);
		quad_on_coarse = nbr_info.quad_on_coarse;
		coarse_id      = nbr_info.id * Side<D>::num_sides + s.opposite().toInt();
		this->own      = true;
		coarse_own     = (nbr_info.ptr != nullptr);
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
		ids.push_back(coarse_id);
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
		idx.push_back(coarse_local_index);
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
		idx.push_back(coarse_global_index);
	}
	void getIfaceTypes(std::deque<IfaceType> &types)
	{
		IfaceType fine_type(IfaceType::fine_to_fine, quad_on_coarse);
		IfaceType coarse_type(IfaceType::fine_to_coarse, quad_on_coarse);
		types.push_back(fine_type);
		types.push_back(coarse_type);
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(nbr_info.rank);
		ranks.push_back(nbr_info.rank);
	}
	void getOwnership(std::deque<bool> &own)
	{
		own.push_back(this->own);
		own.push_back(this->coarse_own);
	}
	void getIncomingProcs(std::set<int> &incoming_procs)
	{
		if (nbr_info.ptr == nullptr) { incoming_procs.insert(nbr_info.rank); }
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index  = rev_map.at(this->id);
		coarse_local_index = rev_map.at(coarse_id);
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index  = rev_map.at(this->local_index);
		coarse_global_index = rev_map.at(coarse_local_index);
	}
};
template <size_t D> class FineIfaceInfo : public IfaceInfo<D>
{
	public:
	FineNbrInfo<D>                                 nbr_info;
	std::array<int, Orthant<D - 1>::num_orthants>  fine_ids;
	std::array<int, Orthant<D - 1>::num_orthants>  fine_local_indexes;
	std::array<int, Orthant<D - 1>::num_orthants>  fine_global_indexes;
	std::array<bool, Orthant<D - 1>::num_orthants> fine_own;
	FineIfaceInfo(Domain<D> &d, Side<D> s)
	{
		nbr_info  = d.getFineNbrInfo(s);
		this->id  = d.id * Side<D>::num_sides + s.toInt();
		this->own = true;
		for (size_t i = 0; i < fine_ids.size(); i++) {
			fine_ids[i] = nbr_info.ids[i] * Side<D>::num_sides + s.opposite().toInt();
			fine_own[i] = (nbr_info.ptrs[i] != nullptr);
		}
	}
	void getIdxAndTypes(std::deque<int> &idx, std::deque<IfaceType> &types)
	{
		idx.push_back(this->local_index);
		types.push_back(IfaceType::coarse_to_coarse);
		for (size_t i = 0; i < fine_local_indexes.size(); i++) {
			idx.push_back(fine_local_indexes[i]);
			IfaceType type(IfaceType::coarse_to_fine, i);
			types.push_back(type);
		}
	}
	void getIds(std::deque<int> &ids)
	{
		ids.push_back(this->id);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			ids.push_back(fine_ids[i]);
		}
	}
	void getLocalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->local_index);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			idx.push_back(fine_local_indexes[i]);
		}
	}
	void getGlobalIndexes(std::deque<int> &idx)
	{
		idx.push_back(this->global_index);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			idx.push_back(fine_global_indexes[i]);
		}
	}
	void getIfaceTypes(std::deque<IfaceType> &types)
	{
		types.push_back(IfaceType::coarse_to_coarse);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			IfaceType type(IfaceType::coarse_to_fine, i);
			types.push_back(type);
		}
	}
	void getRanks(std::deque<int> &ranks)
	{
		ranks.push_back(-1);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			ranks.push_back(nbr_info.ranks[i]);
		}
	}
	void getOwnership(std::deque<bool> &own)
	{
		own.push_back(this->own);
		for (size_t i = 0; i < fine_own.size(); i++) {
			own.push_back(fine_own[i]);
		}
	}
	void getIncomingProcs(std::set<int> &incoming_procs)
	{
		for (size_t i = 0; i < nbr_info.ptrs.size(); i++) {
			if (nbr_info.ptrs[i] == nullptr) { incoming_procs.insert(nbr_info.ranks[i]); }
		}
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		this->local_index = rev_map.at(this->id);
		for (size_t i = 0; i < fine_ids.size(); i++) {
			fine_local_indexes[i] = rev_map.at(fine_ids[i]);
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		this->global_index = rev_map.at(this->local_index);
		for (size_t i = 0; i < fine_local_indexes.size(); i++) {
			fine_global_indexes[i] = rev_map.at(fine_local_indexes[i]);
		}
	}
};
template <size_t D> struct SchurDomain : public Domain<D> {
	std::array<IfaceInfo<D> *, Side<D>::num_sides> iface_info;
	SchurDomain()
	{
		iface_info.fill(nullptr);
	}
	SchurDomain(Domain<D> &d) : Domain<D>(d)
	{
		iface_info.fill(nullptr);

		// create iface objects
		for (Side<D> s : Side<D>::getValues()) {
			if (d.hasNbr(s)) {
				switch (d.getNbrType(s)) {
					case NbrType::Normal:
						getIfaceInfoPtr(s) = new NormalIfaceInfo<D>(d, s);
						break;
					case NbrType::Fine:
						getIfaceInfoPtr(s) = new FineIfaceInfo<D>(d, s);
						break;
					case NbrType::Coarse:
						getIfaceInfoPtr(s) = new CoarseIfaceInfo<D>(d, s);
						break;
				}
			}
		}
	}
	IfaceInfo<D> *&getIfaceInfoPtr(Side<D> s)
	{
		return iface_info[s.toInt()];
	}
	NormalIfaceInfo<D> &getNormalIfaceInfo(Side<D> s)
	{
		return *(NormalIfaceInfo<D> *) iface_info[s.toInt()];
	}
	void enumerateIfaces(std::map<int, IfaceSet<D>> &               ifaces,
	                     std::map<int, std::map<int, IfaceSet<D>>> &off_proc_ifaces,
	                     std::set<int> &                            incoming_procs)
	{
		std::array<int, Side<D>::num_sides> ids;
		for (Side<D> s : Side<D>::getValues()) {
			if (this->hasNbr(s)) {
				ids[s.toInt()] = getIfaceInfoPtr(s)->id;
			} else {
				ids[s.toInt()] = -1;
			}
		}
		std::deque<int>       iface_ids;
		std::deque<IfaceType> iface_types;
		std::deque<Side<D>>   iface_sides;
		std::deque<bool>      iface_own;
		std::deque<int>       iface_ranks;
		for (Side<D> s : Side<D>::getValues()) {
			if (this->hasNbr(s)) {
				int num_added = iface_ids.size();

				getIfaceInfoPtr(s)->getIds(iface_ids);

				num_added = iface_ids.size() - num_added;
				for (int i = 0; i < num_added; i++) {
					iface_sides.push_back(s);
				}

				getIfaceInfoPtr(s)->getIfaceTypes(iface_types);
				getIfaceInfoPtr(s)->getRanks(iface_ranks);
				getIfaceInfoPtr(s)->getOwnership(iface_own);
				getIfaceInfoPtr(s)->getIncomingProcs(incoming_procs);
			}
		}
		for (size_t i = 0; i < iface_ids.size(); i++) {
			int       id   = iface_ids[i];
			IfaceType type = iface_types[i];
			Side<D>   s    = iface_sides[i];
			int       rank = iface_ranks[i];
			if (iface_own[i]) {
				IfaceSet<D> &ifs = ifaces[id];
				ifs.id           = id;
				ifs.insert(Iface<D>(ids, type, s, this->neumann));
			} else {
				IfaceSet<D> &ifs = off_proc_ifaces[rank][id];
				ifs.id           = id;
				ifs.insert(Iface<D>(ids, type, s, this->neumann));
			}
		}
	}
	std::deque<int> getIds()
	{
		std::deque<int> retval;
		for (Side<D> s : Side<D>::getValues()) {
			if (this->hasNbr(s)) { getIfaceInfoPtr(s)->getIds(retval); }
		}
		return retval;
	}
	void setLocalIndexes(const std::map<int, int> &rev_map)
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (this->hasNbr(s)) { getIfaceInfoPtr(s)->setLocalIndexes(rev_map); }
		}
	}
	void setGlobalIndexes(const std::map<int, int> &rev_map)
	{
		for (Side<D> s : Side<D>::getValues()) {
			if (this->hasNbr(s)) { getIfaceInfoPtr(s)->setGlobalIndexes(rev_map); }
		}
	}
	int getIfaceLocalIndex(Side<D> s)
	{
		return iface_info[s.toInt()]->local_index;
	}
};
extern template struct SchurDomain<2>;
extern template struct SchurDomain<3>;
#endif
