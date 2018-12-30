#include "../Utils.h"
#include "catch.hpp"
using namespace std;
TEST_CASE("getSlice getNbrType works", "[Domain]")
{
	vector<double> cube(4 * 4 * 4);
	for (int i = 0; i < 4 * 4 * 4; i++) {
		cube[i] = i;
	}
	{
		auto slice = Utils::getSlice<2>(&cube[0], 4, Side<3>::bottom);
		for (int xi = 0; xi < 4; xi++) {
			for (int yi = 0; yi < 4; yi++) {
				std::array<int, 2> coords = {{xi, yi}};
				CHECK(slice(coords) == xi + yi * 4);
			}
		}
	}
	{
		auto slice = Utils::getSlice<2>(&cube[0], 4, Side<3>::top);
		for (int xi = 0; xi < 4; xi++) {
			for (int yi = 0; yi < 4; yi++) {
				std::array<int, 2> coords = {{xi, yi}};
				CHECK(slice(coords) == 4 * 4 * 3 + xi + yi * 4);
			}
		}
	}
	{
		auto slice = Utils::getSlice<2>(&cube[0], 4, Side<3>::west);
		for (int xi = 0; xi < 4; xi++) {
			for (int yi = 0; yi < 4; yi++) {
				std::array<int, 2> coords = {{xi, yi}};
				CHECK(slice(coords) == xi*4 + yi * 4*4);
			}
		}
	}
	{
		auto slice = Utils::getSlice<2>(&cube[0], 4, Side<3>::east);
		for (int xi = 0; xi < 4; xi++) {
			for (int yi = 0; yi < 4; yi++) {
				std::array<int, 2> coords = {{xi, yi}};
				CHECK(slice(coords) ==3+ xi*4 + yi * 4*4);
			}
		}
	}
	{
		auto slice = Utils::getSlice<2>(&cube[0], 4, Side<3>::south);
		for (int xi = 0; xi < 4; xi++) {
			for (int yi = 0; yi < 4; yi++) {
				std::array<int, 2> coords = {{xi, yi}};
				CHECK(slice(coords) == xi + yi * 4*4);
			}
		}
	}
	{
		auto slice = Utils::getSlice<2>(&cube[0], 4, Side<3>::north);
		for (int xi = 0; xi < 4; xi++) {
			for (int yi = 0; yi < 4; yi++) {
				std::array<int, 2> coords = {{xi, yi}};
				CHECK(slice(coords) == 4*3+xi + yi * 4*4);
			}
		}
	}
}
