#include <Thunderegg/Side.h>
#include "catch.hpp"
using namespace std;
TEST_CASE("Default constructor works", "[Side]")
{
	Side<3> s;
	REQUIRE(s.toInt() == -1);
}
TEST_CASE("int constructor and toInt() works", "[Side]")
{
	{
		Side<3> s(5);
		REQUIRE(s.toInt() == 5);
	}
	{
		Side<3> s(-8);
		REQUIRE(s.toInt() == -8);
	}
}
TEST_CASE("== operator works", "[Side]")
{
	{
		Side<3> s(5);
		Side<3> other(5);
		Side<3> lower(4);
		Side<3> higher(6);
		REQUIRE(!(s == 4));
		REQUIRE(s == 5);
		REQUIRE(!(s == 6));
		REQUIRE(!(s == lower));
		REQUIRE(s == other);
		REQUIRE(!(s == higher));
	}
	{
		Side<3> s(0);
		Side<3> other(0);
		Side<3> lower(-1);
		Side<3> higher(1);
		REQUIRE(!(s == -1));
		REQUIRE(s == 0);
		REQUIRE(!(s == 1));
		REQUIRE(!(s == lower));
		REQUIRE(s == other);
		REQUIRE(!(s == higher));
	}
}
TEST_CASE("< operator works", "[Side]")
{
	{
		Side<3> s(5);
		Side<3> other(5);
		Side<3> lower(4);
		Side<3> higher(6);
		REQUIRE(!(s < lower));
		REQUIRE(!(s < other));
		REQUIRE(s < higher);
	}
	{
		Side<3> s(0);
		Side<3> other(0);
		Side<3> lower(-1);
		Side<3> higher(1);
		REQUIRE(!(s < lower));
		REQUIRE(!(s < other));
		REQUIRE(s < higher);
	}
}
TEST_CASE("Side<3> values are as expected", "[Side]")
{
	{
		Side<3> s(0);
		int     expected = Side<3>::west;
		REQUIRE(s == expected);
	}
	{
		Side<3> s(1);
		int     expected = Side<3>::east;
		REQUIRE(s == expected);
	}
	{
		Side<3> s(2);
		int     expected = Side<3>::south;
		REQUIRE(s == expected);
	}
	{
		Side<3> s(3);
		int     expected = Side<3>::north;
		REQUIRE(s == expected);
	}
	{
		Side<3> s(4);
		int     expected = Side<3>::bottom;
		REQUIRE(s == expected);
	}
	{
		Side<3> s(5);
		int     expected = Side<3>::top;
		REQUIRE(s == expected);
	}
}
TEST_CASE("side.opposite() values are as expected", "[Side]")
{
	{
		Side<3> s        = Side<3>::west;
		Side<3> expected = Side<3>::east;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side<3> s        = Side<3>::east;
		Side<3> expected = Side<3>::west;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side<3> s        = Side<3>::south;
		Side<3> expected = Side<3>::north;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side<3> s        = Side<3>::north;
		Side<3> expected = Side<3>::south;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side<3> s        = Side<3>::bottom;
		Side<3> expected = Side<3>::top;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side<3> s        = Side<3>::top;
		Side<3> expected = Side<3>::bottom;
		REQUIRE(s.opposite() == expected);
	}
}
TEST_CASE("side.isOnLowerAxis() values are as expected", "[Side]")
{
	{
		Side<3> s = Side<3>::west;
		REQUIRE(s.isLowerOnAxis() == true);
	}
	{
		Side<3> s = Side<3>::east;
		REQUIRE(s.isLowerOnAxis() == false);
	}
	{
		Side<3> s = Side<3>::south;
		REQUIRE(s.isLowerOnAxis() == true);
	}
	{
		Side<3> s = Side<3>::north;
		REQUIRE(s.isLowerOnAxis() == false);
	}
	{
		Side<3> s = Side<3>::bottom;
		REQUIRE(s.isLowerOnAxis() == true);
	}
	{
		Side<3> s = Side<3>::top;
		REQUIRE(s.isLowerOnAxis() == false);
	}
}
