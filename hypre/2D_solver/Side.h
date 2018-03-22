#ifndef SIDEENUM_H
#define SIDEENUM_H
enum class Side { north, east, south, west };
inline Side operator++(Side &s,int i)
{
	s = static_cast<Side>((static_cast<int>(s) + 1)%4);
	return s;
}
inline Side operator--(Side &s,int i)
{
	s = static_cast<Side>((static_cast<int>(s) + 3)%4);
	return s;
}
inline Side operator!(const Side s) { return static_cast<Side>((static_cast<int>(s) + 2) % 4); }

enum class Tilt { center, left, right };
#endif
