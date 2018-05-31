#ifndef SIDEENUM_H
#define SIDEENUM_H
#include <array>
#include <iostream>
#include <vector>
enum class Side : char { west, east, south, north, bottom, top };
enum class Oct : char { bsw, bse, bnw, bne, tsw, tse, tnw, tne };
inline std::array<Side, 6> getSideValues(){
    return std::array<Side,6>({{Side::west,Side::east,Side::south,Side::north,Side::bottom,Side::top}});
}
inline std::array<Oct, 4> getOctsOnSide(Side s)
{
	std::array<Oct, 4> retval;
	switch (s) {
		case Side::west:
			retval = {{Oct::bsw, Oct::bnw, Oct::tsw, Oct::tnw}};
			break;
		case Side::east:
			retval = {{Oct::bse, Oct::bne, Oct::tse, Oct::tne}};
			break;
		case Side::south:
			retval = {{Oct::bsw, Oct::bse, Oct::tsw, Oct::tse}};
			break;
		case Side::north:
			retval = {{Oct::bnw, Oct::bne, Oct::tnw, Oct::tne}};
			break;
		case Side::bottom:
			retval = {{Oct::bsw, Oct::bse, Oct::bnw, Oct::bne}};
			break;
		case Side::top:
			retval = {{Oct::tsw, Oct::tse, Oct::tnw, Oct::tne}};
			break;
	}
	return retval;
}
inline bool isENT(Side s) { return static_cast<int>(s)%2==1; }
inline Side operator++(Side &s, int i)
{
	s = static_cast<Side>((static_cast<char>(s) + 1) % 6);
	return s;
}
inline Side &operator++(Side &s)
{
	s = static_cast<Side>(static_cast<char>(s) + 1);
	return s;
}
inline int operator^(const int &a, const Side &b) { return a * 6 + static_cast<int>(b); }
inline int operator%(const Side &a, const int &b) { return static_cast<char>(a) % b; }
inline Side operator~(const Side s)
{
	if (s % 2 == 0) {
		return static_cast<Side>(static_cast<char>(s) + 1);
	} else {
		return static_cast<Side>(static_cast<char>(s) - 1);
	}
}
inline std::ostream &operator<<(std::ostream &os, const Side &s)
{
	switch (s) {
		case Side::west:
			os << "Side::west";
			break;
		case Side::east:
			os << "Side::east";
			break;
		case Side::south:
			os << "Side::south";
			break;
		case Side::north:
			os << "Side::north";
			break;
		case Side::bottom:
			os << "Side::bottom";
			break;
		case Side::top:
			os << "Side::top";
			break;
	}
	return os;
}
#endif
