#ifndef DOMAIN_H
#define DOMAIN_H
#include "Side.h"
#include <array>
#include <bitset>
#include <map>
#include <vector>
enum class NbrType { Normal, Coarse, Fine };
class NbrInfo
{
	public:
	virtual NbrType getNbrType()                      = 0;
	virtual void getNbrIds(std::vector<int> &nbr_ids) = 0;
	virtual void setGlobalIndexes(std::map<int, int> &rev_map) = 0;
	virtual void setLocalIndexes(std::map<int, int> &rev_map)  = 0;
};
class NormalNbrInfo : public NbrInfo
{
	public:
	int id           = 0;
	int local_index  = 0;
	int global_index = 0;
	NormalNbrInfo(int id) { this->id = id; }
	NbrType           getNbrType() { return NbrType::Normal; }
	void getNbrIds(std::vector<int> &nbr_ids) { nbr_ids.push_back(id); };
	void setGlobalIndexes(std::map<int, int> &rev_map) { global_index = rev_map.at(local_index); }
	void setLocalIndexes(std::map<int, int> &rev_map) { local_index = rev_map.at(id); }
};
class CoarseNbrInfo : public NbrInfo
{
	public:
	int id;
	int local_index;
	int global_index;
	int quad_on_coarse;
	CoarseNbrInfo(int id, int quad_on_coarse)
	{
		this->id             = id;
		this->quad_on_coarse = quad_on_coarse;
	}
	NbrType getNbrType() { return NbrType::Coarse; }
	void getNbrIds(std::vector<int> &nbr_ids) { nbr_ids.push_back(id); };
	void setGlobalIndexes(std::map<int, int> &rev_map) { global_index = rev_map.at(local_index); }
	void setLocalIndexes(std::map<int, int> &rev_map) { local_index = rev_map.at(id); }
};
class FineNbrInfo : public NbrInfo
{
	public:
	std::array<int, 4> ids;
	std::array<int, 4> global_indexes;
	std::array<int, 4> local_indexes;
	FineNbrInfo(std::array<int, 4> ids) { this->ids = ids; }
	NbrType getNbrType() { return NbrType::Fine; }
	void getNbrIds(std::vector<int> &nbr_ids)
	{
		for (int i = 0; i < 4; i++) {
			nbr_ids.push_back(ids[i]);
		}
	};
	void setGlobalIndexes(std::map<int, int> &rev_map)
	{
		for (int i = 0; i < 4; i++) {
			global_indexes[i] = rev_map.at(local_indexes[i]);
		}
	}
	void setLocalIndexes(std::map<int, int> &rev_map)
	{
		for (int i = 0; i < 4; i++) {
			local_indexes[i] = rev_map.at(ids[i]);
		}
	}
};
/**
 * @brief A structure that represents a domain and its relation to other domains.
 */
struct Domain {
	/**
	 * @brief The domain's own id
	 */
	int id        = 0;
	int id_local  = 0;
	int id_global = 0;
	int n         = 10;

	int refine_level = 1;
	int parent_id    = -1;

	std::array<NbrInfo *, 6> nbr_info = {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}};
	std::array<int, 8>       child_id = {{-1, -1, -1, -1, -1, -1, -1, -1}};

	std::bitset<6> neumann;
	bool           zero_patch = false;
	/**
	 * @brief The lower left x coordinate of domain
	 */
	double x_start = 0;
	/**
	 * @brief The lower left y coordinate of domain
	 */
	double y_start = 0;
	double z_start = 0;
	/**
	 * @brief length of domain in x direction
	 */
	double x_length = 1;
	/**
	 * @brief length of odmain in y direction
	 */
	double y_length = 1;
	double z_length = 1;

	friend bool operator<(const Domain &l, const Domain &r) { return l.id < r.id; }
	NbrInfo *&getNbrInfoPtr(Side s) { return nbr_info[static_cast<int>(s)]; }
	NbrType getNbrType(Side s) { return nbr_info[static_cast<int>(s)]->getNbrType(); }
	NormalNbrInfo &getNormalNbrInfo(Side s)
	{
		return *(NormalNbrInfo *) nbr_info[static_cast<int>(s)];
	}
	CoarseNbrInfo &getCoarseNbrInfo(Side s)
	{
		return *(CoarseNbrInfo *) nbr_info[static_cast<int>(s)];
	}
	FineNbrInfo &getFineNbrInfo(Side s) { return *(FineNbrInfo *) nbr_info[static_cast<int>(s)]; }
	inline bool hasNbr(Side s) const { return nbr_info[static_cast<int>(s)] != nullptr; }
	inline bool isNeumann(Side s) const { return neumann[static_cast<int>(s)]; }
	void setLocalNeighborIndexes(std::map<int, int> &rev_map)
	{
		id_local = rev_map.at(id);
		for (Side s : getSideValues()) {
			if (hasNbr(s)) {
				getNbrInfoPtr(s)->setLocalIndexes(rev_map);
			}
		}
	}

	void setGlobalNeighborIndexes(std::map<int, int> &rev_map)
	{
		id_global = rev_map.at(id_local);
		for (Side s : getSideValues()) {
			if (hasNbr(s)) {
				getNbrInfoPtr(s)->setGlobalIndexes(rev_map);
			}
		}
	}
	void setNeumann()
	{
		for (int q = 0; q < 6; q++) {
			neumann[q] = !hasNbr(static_cast<Side>(q));
		}
	}
	std::vector<int> getNbrIds()
	{
		std::vector<int> retval;
		for (Side s : getSideValues()) {
			if (hasNbr(s)) {
				getNbrInfoPtr(s)->getNbrIds(retval);
			}
		}
		return retval;
	}
};
#endif
