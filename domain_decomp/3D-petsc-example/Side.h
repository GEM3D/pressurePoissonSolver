#ifndef SIDEENUM_H
#define SIDEENUM_H
#include <iostream>
enum class Side : char { west, east, south, north, bottom, top };
inline Side operator++(Side &s, int i)
{
	s = static_cast<Side>((static_cast<char>(s) + 1) % 6);
	return s;
}
inline int operator^(const int &a, const Side &b) { return a * 6 + static_cast<char>(b); }
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
