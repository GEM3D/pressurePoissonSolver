#ifndef DOMAIN_H
#define DOMAIN_H
#include "BufferWriter.h"
#include "Side.h"
#include <array>
#include <bitset>
#include <map>
#include <vector>
enum class NbrType { Normal, Coarse, Fine };
struct NbrInfo;
struct NormalNbrInfo;
struct CoarseNbrInfo;
struct FineNbrInfo;
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

	std::array<int, 8> child_id = {{-1, -1, -1, -1, -1, -1, -1, -1}};

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
	 * @brief length of domain in y direction
	 */
	double                   y_length = 1;
	double                   z_length = 1;
	std::array<NbrInfo *, 6> nbr_info = {{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr}};

	~Domain();
	friend bool      operator<(const Domain &l, const Domain &r) { return l.id < r.id; }
	NbrInfo *&       getNbrInfoPtr(Side s);
	NbrType          getNbrType(Side s);
	NormalNbrInfo &  getNormalNbrInfo(Side s);
	CoarseNbrInfo &  getCoarseNbrInfo(Side s);
	FineNbrInfo &    getFineNbrInfo(Side s);
	inline bool      hasNbr(Side s) const;
	inline bool      isNeumann(Side s) const;
	void             setLocalNeighborIndexes(std::map<int, int> &rev_map);
	void             setGlobalNeighborIndexes(std::map<int, int> &rev_map);
	void             setNeumann();
	std::vector<int> getNbrIds();
	int              serialize(char *buffer);
	int              serializeSize();
	static Domain    deserialize(char *buffer);
	void             setPtrs(std::map<int, Domain> &domains);
	void             updateRank(int rank);
};
class NbrInfo
{
	public:
	virtual ~NbrInfo()                                            = default;
	virtual NbrType getNbrType()                                  = 0;
	virtual void    getNbrIds(std::vector<int> &nbr_ids)          = 0;
	virtual void    setGlobalIndexes(std::map<int, int> &rev_map) = 0;
	virtual void    setLocalIndexes(std::map<int, int> &rev_map)  = 0;
	virtual void    setPtrs(std::map<int, Domain> &domains)       = 0;
	virtual void    updateRank(int new_rank)                      = 0;
	virtual void    updateRankOnNeighbors(int new_rank, Side s)   = 0;
};
class NormalNbrInfo : public NbrInfo
{
	public:
	Domain *ptr          = nullptr;
	int     rank         = 0;
	int     id           = 0;
	int     local_index  = 0;
	int     global_index = 0;
	NormalNbrInfo() {}
	~NormalNbrInfo() = default;
	NormalNbrInfo(int id) { this->id = id; }
	NbrType getNbrType() { return NbrType::Normal; }
	void    getNbrIds(std::vector<int> &nbr_ids) { nbr_ids.push_back(id); };
	void setGlobalIndexes(std::map<int, int> &rev_map) { global_index = rev_map.at(local_index); }
	void setLocalIndexes(std::map<int, int> &rev_map) { local_index = rev_map.at(id); }
	void setPtrs(std::map<int, Domain> &domains)
	{
		try {
			ptr = &domains.at(id);
		} catch (std::out_of_range) {
			ptr = nullptr;
		}
	}
	void updateRank(int new_rank) { rank = new_rank; }
	void updateRankOnNeighbors(int new_rank, Side s)
	{
		ptr->getNbrInfoPtr(~s)->updateRank(new_rank);
	}
};
class CoarseNbrInfo : public NbrInfo
{
	public:
	Domain *ptr  = nullptr;
	int     rank = 0;
	int     id;
	int     local_index;
	int     global_index;
	int     quad_on_coarse;
	CoarseNbrInfo() {}
	~CoarseNbrInfo() = default;
	CoarseNbrInfo(int id, int quad_on_coarse)
	{
		this->id             = id;
		this->quad_on_coarse = quad_on_coarse;
	}
	NbrType getNbrType() { return NbrType::Coarse; }
	void    getNbrIds(std::vector<int> &nbr_ids) { nbr_ids.push_back(id); };
	void setGlobalIndexes(std::map<int, int> &rev_map) { global_index = rev_map.at(local_index); }
	void setLocalIndexes(std::map<int, int> &rev_map) { local_index = rev_map.at(id); }
	void setPtrs(std::map<int, Domain> &domains)
	{
		try {
			ptr = &domains.at(id);
		} catch (std::out_of_range) {
			ptr = nullptr;
		}
	}
	void updateRank(int new_rank) { rank = new_rank; }
	void updateRankOnNeighbors(int new_rank, Side s)
	{
		ptr->getNbrInfoPtr(~s)->updateRank(new_rank);
	}
};
class FineNbrInfo : public NbrInfo
{
	public:
	std::array<Domain *, 4> ptrs  = {nullptr, nullptr, nullptr, nullptr};
	std::array<int, 4>      ranks = {0, 0, 0, 0};
	std::array<int, 4>      ids;
	std::array<int, 4>      global_indexes;
	std::array<int, 4>      local_indexes;
	FineNbrInfo() {}
	~FineNbrInfo() = default;
	FineNbrInfo(std::array<int, 4> ids) { this->ids = ids; }
	NbrType getNbrType() { return NbrType::Fine; }
	void    getNbrIds(std::vector<int> &nbr_ids)
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
	void setPtrs(std::map<int, Domain> &domains)
	{
		for (int i = 0; i < 4; i++) {
			try {
				ptrs[i] = &domains[ids[i]];
			} catch (std::out_of_range) {
				ptrs[i] = nullptr;
			}
		}
	}
	void updateRank(int new_rank) { ranks.fill(new_rank); }

	void updateRankOnNeighbors(int new_rank, Side s)
	{
		for (int i = 0; i < 4; i++) {
			ptrs[i]->getNbrInfoPtr(~s)->updateRank(new_rank);
		}
	}
};
inline Domain::~Domain()
{
    /*
	for (NbrInfo *info : nbr_info) {
		delete info;
	}
    */
}
inline NbrInfo *&Domain::getNbrInfoPtr(Side s) { return nbr_info[static_cast<int>(s)]; }
inline NbrType   Domain::getNbrType(Side s) { return nbr_info[static_cast<int>(s)]->getNbrType(); }
inline NormalNbrInfo &Domain::getNormalNbrInfo(Side s)
{
	return *(NormalNbrInfo *) nbr_info[static_cast<int>(s)];
}
inline CoarseNbrInfo &Domain::getCoarseNbrInfo(Side s)
{
	return *(CoarseNbrInfo *) nbr_info[static_cast<int>(s)];
}
inline FineNbrInfo &Domain::getFineNbrInfo(Side s)
{
	return *(FineNbrInfo *) nbr_info[static_cast<int>(s)];
}
inline bool Domain::hasNbr(Side s) const { return nbr_info[static_cast<int>(s)] != nullptr; }
inline bool Domain::isNeumann(Side s) const { return neumann[static_cast<int>(s)]; }
inline void Domain::setLocalNeighborIndexes(std::map<int, int> &rev_map)
{
	id_local = rev_map.at(id);
	for (Side s : getSideValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setLocalIndexes(rev_map); }
	}
}

inline void Domain::setGlobalNeighborIndexes(std::map<int, int> &rev_map)
{
	id_global = rev_map.at(id_local);
	for (Side s : getSideValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setGlobalIndexes(rev_map); }
	}
}
inline void Domain::setNeumann()
{
	for (int q = 0; q < 6; q++) {
		neumann[q] = !hasNbr(static_cast<Side>(q));
	}
}
inline std::vector<int> Domain::getNbrIds()
{
	std::vector<int> retval;
	for (Side s : getSideValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->getNbrIds(retval); }
	}
	return retval;
}

inline int Domain::serialize(char *buffer)
{
	BufferWriter writer(buffer);
	writer << id;
	writer << n;
	writer << refine_level;
	writer << parent_id;
	writer << child_id;
	writer << neumann;
	writer << zero_patch;
	writer << x_start;
	writer << y_start;
	writer << z_start;
	writer << x_length;
	writer << y_length;
	writer << z_length;
	std::bitset<6> has_nbr;
	for (int i = 0; i < 6; i++) {
		has_nbr[i] = nbr_info[i] != nullptr;
	}
	writer << has_nbr;
	for (Side s : getSideValues()) {
		if (hasNbr(s)) {
			NbrType type = getNbrType(s);
			writer << type;
			switch (type) {
				case NbrType::Normal:
					writer << getNormalNbrInfo(s);
					break;
				case NbrType::Fine:
					writer << getFineNbrInfo(s);
					break;
				case NbrType::Coarse:
					writer << getCoarseNbrInfo(s);
					break;
			}
		}
	}
	return writer.getPos();
}
inline Domain Domain::deserialize(char *buffer)
{
	BufferReader reader(buffer);
	Domain       d;
	reader >> d.id;
	reader >> d.n;
	reader >> d.refine_level;
	reader >> d.parent_id;
	reader >> d.child_id;
	reader >> d.neumann;
	reader >> d.zero_patch;
	reader >> d.x_start;
	reader >> d.y_start;
	reader >> d.z_start;
	reader >> d.x_length;
	reader >> d.y_length;
	reader >> d.z_length;
	std::bitset<6> has_nbr;
	reader >> has_nbr;
	for (int i = 0; i < 6; i++) {
		if (has_nbr[i]) {
			NbrType type;
			reader >> type;
			NbrInfo *info = nullptr;
			switch (type) {
				case NbrType::Normal:
					info = new NormalNbrInfo();
					reader >> *(NormalNbrInfo *) info;
					break;
				case NbrType::Fine:
					info = new FineNbrInfo();
					reader >> *(FineNbrInfo *) info;
					break;
				case NbrType::Coarse:
					info = new CoarseNbrInfo();
					reader >> *(CoarseNbrInfo *) info;
					break;
			}
			d.nbr_info[i] = info;
		}
	}
	return d;
}
inline void Domain::setPtrs(std::map<int, Domain> &domains)
{
	for (Side s : getSideValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->setPtrs(domains); }
	}
}
inline void Domain::updateRank(int rank)
{
	for (Side s : getSideValues()) {
		if (hasNbr(s)) { getNbrInfoPtr(s)->updateRankOnNeighbors(rank, s); }
	}
}
#endif
