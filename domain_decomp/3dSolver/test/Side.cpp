#include "Side.h"
#include "catch.hpp"
using namespace std;
TEST_CASE("Default constructor works", "[Side]")
{
	Side s;
	REQUIRE(s.toInt() == -1);
}
TEST_CASE("int constructor and toInt() works", "[Side]")
{
	{
		Side s(5);
		REQUIRE(s.toInt() == 5);
	}
	{
		Side s(-8);
		REQUIRE(s.toInt() == -8);
	}
}
TEST_CASE("== operator works", "[Side]")
{
	{
		Side s(5);
		Side other(5);
		Side lower(4);
		Side higher(6);
		REQUIRE(!(s == 4));
		REQUIRE(s == 5);
		REQUIRE(!(s == 6));
		REQUIRE(!(s == lower));
		REQUIRE(s == other);
		REQUIRE(!(s == higher));
	}
	{
		Side s(0);
		Side other(0);
		Side lower(-1);
		Side higher(1);
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
		Side s(5);
		Side other(5);
		Side lower(4);
		Side higher(6);
		REQUIRE(!(s < lower));
		REQUIRE(!(s < other));
		REQUIRE(s < higher);
	}
	{
		Side s(0);
		Side other(0);
		Side lower(-1);
		Side higher(1);
		REQUIRE(!(s < lower));
		REQUIRE(!(s < other));
		REQUIRE(s < higher);
	}
}
TEST_CASE("Side values are as expected", "[Side]")
{
	{
		Side s(0);
		int  expected = Side::west;
		REQUIRE(s == expected);
	}
	{
		Side s(1);
		int  expected = Side::east;
		REQUIRE(s == expected);
	}
	{
		Side s(2);
		int  expected = Side::south;
		REQUIRE(s == expected);
	}
	{
		Side s(3);
		int  expected = Side::north;
		REQUIRE(s == expected);
	}
	{
		Side s(4);
		int  expected = Side::bottom;
		REQUIRE(s == expected);
	}
	{
		Side s(5);
		int  expected = Side::top;
		REQUIRE(s == expected);
	}
}
TEST_CASE("side.opposite() values are as expected", "[Side]")
{
	{
		Side s        = Side::west;
		Side expected = Side::east;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side s        = Side::east;
		Side expected = Side::west;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side s        = Side::south;
		Side expected = Side::north;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side s        = Side::north;
		Side expected = Side::south;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side s        = Side::bottom;
		Side expected = Side::top;
		REQUIRE(s.opposite() == expected);
	}
	{
		Side s        = Side::top;
		Side expected = Side::bottom;
		REQUIRE(s.opposite() == expected);
	}
}
TEST_CASE("side.isOnLowerAxis() values are as expected", "[Side]")
{
	{
		Side s = Side::west;
		REQUIRE(s.isLowerOnAxis() == true);
	}
	{
		Side s = Side::east;
		REQUIRE(s.isLowerOnAxis() == false);
	}
	{
		Side s = Side::south;
		REQUIRE(s.isLowerOnAxis() == true);
	}
	{
		Side s = Side::north;
		REQUIRE(s.isLowerOnAxis() == false);
	}
	{
		Side s = Side::bottom;
		REQUIRE(s.isLowerOnAxis() == true);
	}
	{
		Side s = Side::top;
		REQUIRE(s.isLowerOnAxis() == false);
	}
}
