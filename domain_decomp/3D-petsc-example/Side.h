#ifndef SIDEENUM_H
#define SIDEENUM_H
enum class Side { west, east, south, north, bottom, top };
inline Side operator++(Side &s, int i)
{
	s = static_cast<Side>((static_cast<int>(s) + 1) % 6);
	return s;
}
inline Side operator+(const Side &a, const Side &b)
{
	Side s;
	s = static_cast<Side>((static_cast<int>(a) + static_cast<int>(b)) % 4);
	return s;
}
inline Side operator-(const Side &a, const Side &b)
{
	Side s;
	s = static_cast<Side>((static_cast<int>(a) - static_cast<int>(b)) % 4);
	return s;
}
inline Side operator+(const Side &a, const int &b)
{
	Side s;
	s = static_cast<Side>((static_cast<int>(a) + b) % 4);
	return s;
}
inline int operator+(const int &a, const Side &b) { return a + static_cast<int>(b); }
inline Side operator--(Side &s, int i)
{
	s = static_cast<Side>((static_cast<int>(s) + 3) % 4);
	return s;
}
inline int operator^(const int &a, const Side &b) { return a * 6 + static_cast<int>(b); }
inline int operator%(const Side &a, const int &b) { return static_cast<int>(a) % b; }
inline Side operator~(const Side s)
{
	if (s % 2 == 0) {
		return static_cast<Side>(static_cast<int>(s) + 1);
	} else {
		return static_cast<Side>(static_cast<int>(s) - 1);
	}
}
enum class Tilt { center, left, right };
#endif
