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
/**
 * \file
 */

#ifndef THUNDEREGG_PATCHINFO_H
#define THUNDEREGG_PATCHINFO_H
#include <Thunderegg/BufferWriter.h>
#include <Thunderegg/Serializable.h>
#include <Thunderegg/Side.h>
#include <Thunderegg/TypeDefs.h>
#include <array>
#include <bitset>
#include <map>
#include <memory>
#include <vector>

/**
 * @brief The type of neighbor
 */
enum class NbrType {
	/**
	 * @brief The neighbor is at the same refinement level.
	 */
	Normal,
	/**
	 * @brief The neighbor is at a coarser refinement level.
	 */
	Coarse,
	/**
	 * @brief The nighbor is at a finer refinement level.
	 */
	Fine
};

template <size_t D> struct NbrInfo;
template <size_t D> struct NormalNbrInfo;
template <size_t D> struct CoarseNbrInfo;
template <size_t D> struct FineNbrInfo;

/**
 * @brief Contains useful information for a patch
 *
 * This contains information for a specific patch. Information like:
 * * The globally unique id of this patch
 * * The local and global indexes of this patch in the domain
 * * The parent patch (if there is one)
 *
 * It also contains information for a patch's neighbor:
 * * What are the neighbors id?
 * * Are the neighbors at the same refinement level? Are coarser or finer?
 *
 * @tparam D the number of cartesian dimensions in the patch
 */
template <size_t D> struct PatchInfo : public Serializable {
	/**
	 * @brief The globally unique ID of the patch
	 * This ID only needs to be unique within a Domain.
	 */
	int id = 0;
	/**
	 * @brief The local index of the patch in the Domain.
	 */
	int local_index = 0;
	/**
	 * @brief The global index of the patch in the Domain.
	 */
	int global_index = 0;
	/**
	 * @brief The refinement level
	 */
	int refine_level = 1;
	/**
	 * @brief The id of the parent patch.
	 *
	 * Set to -1 if there is no parent.
	 */
	int parent_id = -1;
	/**
	 * @brief The orthant of the parent that this parent resides on.
	 */
	Orthant<D> orth_on_parent;
	/**
	 * @brief Whether the patch has neumann boundary conditions on one side.
	 */
	std::bitset<Side<D>::num_sides> neumann;
	/**
	 * @brief The number of cells in each direction
	 */
	std::array<int, D> ns;
	/**
	 * @brief The lower-left-bottom index of the patch
	 */
	std::array<double, D> starts;
	/**
	 * @brief The cell spacings in each direction
	 */
	std::array<double, D> spacings;
	/**
	 * @brief Nbr info objects for each side.
	 * If there is no neighbor, it shoud be set to nullptr.
	 */
	std::array<NbrInfo<D> *, Side<D>::num_sides> nbr_info;

	/**
	 * @brief Construct a new Patch Info object
	 * starts, ns, and spacings are all set to 0
	 */
	PatchInfo()
	{
		starts.fill(0);
		nbr_info.fill(nullptr);
		ns.fill(0);
		spacings.fill(0);
	}
	/**
	 * @brief Destroy the Patch Info object
	 */
	~PatchInfo();
	/**
	 * @brief Compare the ids of the patches
	 *
	 * @param l left operand
	 * @param r right operand
	 * @return true if r's id is lower
	 * @return false if r's id is not lower
	 */
	friend bool operator<(const PatchInfo &l, const PatchInfo &r)
	{
		return l.id < r.id;
	}
	/**
	 * @brief Get the NbrInfo pointer for a given side
	 *
	 * @param s the side
	 * @return NbrInfo<D>*& the NbrInfo ptr
	 */
	NbrInfo<D> *&getNbrInfoPtr(Side<D> s);
	/**
	 * @brief Get the NbrType for a side
	 *
	 * @param s the side
	 * @return The NbrType
	 */
	NbrType getNbrType(Side<D> s) const;
	/**
	 * @brief Get the NormalNbrInfo object for a side
	 *
	 * Neighbor must be of Normal type, otherwise behavior is undefined.
	 *
	 * @param s the side
	 * @return NormalNbrInfo<D>& the object
	 */
	NormalNbrInfo<D> &getNormalNbrInfo(Side<D> s) const;
	/**
	 * @brief Get the CoarseNbrInfo object
	 *

	 * @param s the side
	 * @return CoarseNbrInfo<D>& the object
	*/
	CoarseNbrInfo<D> &getCoarseNbrInfo(Side<D> s) const;
	/**
	 * @brief Get the FineNbrInfo object
	 *
	 * Neighbor must be of Fine type, otherwise behavior is undefined.
	 *
	 * @param s the side
	 * @return FineNbrInfo<D>& the object
	 */
	FineNbrInfo<D> &getFineNbrInfo(Side<D> s) const;
	/**
	 * @brief Return whether the patch has a neighbor
	 *
	 * @param s the side
	 * @return true if the is neighbor
	 * @return false if at domain boundary
	 */
	inline bool hasNbr(Side<D> s) const;
	/**
	 * @brief Return whether the patch has a coarser parent
	 *
	 * @return true if there is a parent
	 */
	inline bool hasCoarseParent() const;
	/**
	 * @brief Return wheter the boundary conditions are neumann
	 *
	 * @param s the side
	 */
	inline bool isNeumann(Side<D> s) const;
	/**
	 * @brief Set the local indexes in the NbrInfo objects
	 *
	 * @param rev_map map from id to local_index
	 */
	void setLocalNeighborIndexes(std::map<int, int> &rev_map);
	/**
	 * @brief Set the global indexes in the NbrInfo objects
	 *
	 * @param rev_map map form local_index to global_index
	 */
	void setGlobalNeighborIndexes(std::map<int, int> &rev_map);
	/**
	 * @brief Set the neumann boundary conditions
	 *
	 * @param inf the function for determining boudnary conditions
	 */
	void setNeumann(IsNeumannFunc<D> inf);
	/**
	 * @brief return a vector of neighbor ids
	 */
	std::vector<int> getNbrIds();
	/**
	 * @brief set the ptrs to PatchInfo objects in the NbrInfo objects.
	 *
	 * @param patches map of id to patches on this processor
	 */
	void setPtrs(std::map<int, std::shared_ptr<PatchInfo>> &patches);
	/**
	 * @brief change the rank of this object and update the neighbors. (if they are on the same
	 * processor)
	 *
	 * @param rank the mpi rank
	 */
	void updateRank(int rank);
	int  serialize(char *buffer) const;
	int  deserialize(char *buffer);
};
template <size_t D> class NbrInfo : virtual public Serializable
{
	public:
	virtual ~NbrInfo()                                                             = default;
	virtual NbrType getNbrType()                                                   = 0;
	virtual void    getNbrIds(std::vector<int> &nbr_ids)                           = 0;
	virtual void    setGlobalIndexes(std::map<int, int> &rev_map)                  = 0;
	virtual void    setLocalIndexes(std::map<int, int> &rev_map)                   = 0;
	virtual void    setPtrs(std::map<int, std::shared_ptr<PatchInfo<D>>> &domains) = 0;
	virtual void    updateRankOnNeighbors(int new_rank, Side<D> s)                 = 0;
};
template <size_t D> class NormalNbrInfo : public NbrInfo<D>
{
	public:
	std::shared_ptr<PatchInfo<D>> ptr          = nullptr;
	int                           rank         = 0;
	int                           id           = 0;
	int                           local_index  = 0;
	int                           global_index = 0;
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
	void setPtrs(std::map<int, std::shared_ptr<PatchInfo<D>>> &domains)
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
	std::shared_ptr<PatchInfo<D>> ptr;
	int                           rank = 0;
	int                           id;
	int                           local_index;
	int                           global_index;
	int                           quad_on_coarse;
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
	void setPtrs(std::map<int, std::shared_ptr<PatchInfo<D>>> &domains)
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
	std::array<std::shared_ptr<PatchInfo<D>>, Orthant<D - 1>::num_orthants> ptrs;
	std::array<int, Orthant<D - 1>::num_orthants>                           ranks;
	std::array<int, Orthant<D - 1>::num_orthants>                           ids;
	std::array<int, Orthant<D - 1>::num_orthants>                           global_indexes;
	std::array<int, Orthant<D - 1>::num_orthants>                           local_indexes;
	FineNbrInfo()
	{
		ptrs.fill(nullptr);
		ranks.fill(0);
	}
	~FineNbrInfo() = default;
	FineNbrInfo(std::array<int, Orthant<D>::num_orthants / 2> ids)
	{
		ptrs.fill(nullptr);
		ranks.fill(0);
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
	void setPtrs(std::map<int, std::shared_ptr<PatchInfo<D>>> &domains)
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
template <size_t D> inline PatchInfo<D>::~PatchInfo()
{
	/*
	for (NbrInfo *info : nbr_info) {
	    delete info;
	}
	*/
}
template <size_t D> inline NbrInfo<D> *&PatchInfo<D>::getNbrInfoPtr(Side<D> s)
{
	return nbr_info[s.toInt()];
}
template <size_t D> inline NbrType PatchInfo<D>::getNbrType(Side<D> s) const
{
	return nbr_info[s.toInt()]->getNbrType();
}
template <size_t D> inline NormalNbrInfo<D> &PatchInfo<D>::getNormalNbrInfo(Side<D> s) const
{
	return *(NormalNbrInfo<D> *) nbr_info[s.toInt()];
}
template <size_t D> inline CoarseNbrInfo<D> &PatchInfo<D>::getCoarseNbrInfo(Side<D> s) const
{
	return *(CoarseNbrInfo<D> *) nbr_info[s.toInt()];
}
template <size_t D> inline FineNbrInfo<D> &PatchInfo<D>::getFineNbrInfo(Side<D> s) const
{
	return *(FineNbrInfo<D> *) nbr_info[s.toInt()];
}
template <size_t D> inline bool PatchInfo<D>::hasNbr(Side<D> s) const
{
	return nbr_info[s.toInt()] != nullptr;
}
template <size_t D> inline bool PatchInfo<D>::hasCoarseParent() const
{
	return orth_on_parent.toInt() != -1;
}
template <size_t D> inline bool PatchInfo<D>::isNeumann(Side<D> s) const
{
	return neumann[s.toInt()];
}
template <size_t D> inline void PatchInfo<D>::setLocalNeighborIndexes(std::map<int, int> &rev_map)
{
	local_index = rev_map.at(id);
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setLocalIndexes(rev_map); }
	}
}
template <size_t D> inline void PatchInfo<D>::setGlobalNeighborIndexes(std::map<int, int> &rev_map)
{
	global_index = rev_map.at(local_index);
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setGlobalIndexes(rev_map); }
	}
}
template <size_t D> inline void PatchInfo<D>::setNeumann(IsNeumannFunc<D> inf)
{
	for (Side<D> s : Side<D>::getValues()) {
		if (!hasNbr(s)) {
			std::array<double, D> bound_start = starts;
			if (!s.isLowerOnAxis()) { bound_start[s.axis()] += spacings[s.axis()] * ns[s.axis()]; }
			std::array<double, D> bound_end = bound_start;
			for (size_t dir = 0; dir < D; dir++) {
				if (dir != s.axis()) { bound_end[dir] += spacings[dir] * ns[dir]; }
			}
			neumann[s.toInt()] = inf(s, bound_end, bound_start);
		}
	}
}
template <size_t D> inline std::vector<int> PatchInfo<D>::getNbrIds()
{
	std::vector<int> retval;
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->getNbrIds(retval); }
	}
	return retval;
}

template <size_t D> inline int PatchInfo<D>::serialize(char *buffer) const
{
	BufferWriter writer(buffer);
	writer << id;
	writer << ns;
	writer << refine_level;
	writer << parent_id;
	writer << orth_on_parent;
	writer << neumann;
	writer << starts;
	writer << spacings;
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
template <size_t D> inline int PatchInfo<D>::deserialize(char *buffer)
{
	BufferReader reader(buffer);
	reader >> id;
	reader >> ns;
	reader >> refine_level;
	reader >> parent_id;
	reader >> orth_on_parent;
	reader >> neumann;
	reader >> starts;
	reader >> spacings;
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
template <size_t D>
inline void PatchInfo<D>::setPtrs(std::map<int, std::shared_ptr<PatchInfo>> &domains)
{
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setPtrs(domains); }
	}
}
template <size_t D> inline void PatchInfo<D>::updateRank(int rank)
{
	rank = rank;
	for (Side<D> s : Side<D>::getValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->updateRankOnNeighbors(rank, s); }
	}
}
extern template class PatchInfo<2>;
extern template class PatchInfo<3>;
#endif
