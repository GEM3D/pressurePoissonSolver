#ifndef INTERPCASE_H
#define INTERPCASE_H
enum class InterpCase {
	normal,
	coarse_from_coarse,
	coarse_from_fine_on_left,
	coarse_from_fine_on_right,
	fine_from_fine_on_left,
	fine_from_fine_on_right,
	fine_from_coarse_to_fine_on_left,
	fine_from_coarse_to_fine_on_right
};
#define INTERPCASE_H
#endif
