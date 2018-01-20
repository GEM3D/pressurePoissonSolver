#ifndef IFACETYPE_H
#define IFACETYPE_H
enum class IfaceType {
	normal,
	fine_to_fine_0,
	fine_to_fine_1,
	fine_to_fine_2,
	fine_to_fine_3,
	fine_to_coarse_0,
	fine_to_coarse_1,
	fine_to_coarse_2,
	fine_to_coarse_3,
	coarse_to_coarse,
	coarse_to_fine_0,
	coarse_to_fine_1,
	coarse_to_fine_2,
	coarse_to_fine_3
};
inline IfaceType operator+(const IfaceType &a, const int &b)
{
	return static_cast<IfaceType>(static_cast<int>(a) + b);
}
#endif
