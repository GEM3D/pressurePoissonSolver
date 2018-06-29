#include "Side.h"
#include "catch.hpp"
using namespace std;
TEST_CASE("Octant Default constructor works", "[Octant]")
{
	Octant o;
	REQUIRE(o.toInt() == -1);
}
TEST_CASE("Octant int constructor and toInt() works", "[Octant]")
{
	{
		Octant o(5);
		REQUIRE(o.toInt() == 5);
	}
	{
		Octant o(-8);
		REQUIRE(o.toInt() == -8);
	}
}
TEST_CASE("Otant == operator works", "[Octant]")
{
	{
		Octant o(5);
		REQUIRE(!(o == 4));
		REQUIRE(o == 5);
		REQUIRE(!(o == 6));
	}
	{
		Side o(0);
		REQUIRE(!(o == -1));
		REQUIRE(o == 0);
		REQUIRE(!(o == 1));
	}
}
TEST_CASE("Octant values are as expected", "[Octant]")
{
	{
		Octant o(0);
		Octant expected = Octant::bsw;
		REQUIRE(o == expected);
	}
	{
		Octant o(1);
		Octant expected = Octant::bse;
		REQUIRE(o == expected);
	}
	{
		Octant o(2);
		Octant expected = Octant::bnw;
		REQUIRE(o == expected);
	}
	{
		Octant o(3);
		Octant expected = Octant::bne;
		REQUIRE(o == expected);
	}
	{
		Octant o(4);
		Octant expected = Octant::tsw;
		REQUIRE(o == expected);
	}
	{
		Octant o(5);
		Octant expected = Octant::tse;
		REQUIRE(o == expected);
	}
	{
		Octant o(6);
		Octant expected = Octant::tnw;
		REQUIRE(o == expected);
	}
	{
		Octant o(7);
		Octant expected = Octant::tne;
		REQUIRE(o == expected);
	}
}
TEST_CASE("Octant getInteriorNbrOnSide() is as expected", "[Octant]")
{
	{
		Octant o = Octant::bsw;
		REQUIRE(o.getInteriorNbrOnSide(Side::east) == Octant(Octant::bse));
		REQUIRE(o.getInteriorNbrOnSide(Side::north) == Octant(Octant::bnw));
		REQUIRE(o.getInteriorNbrOnSide(Side::top) == Octant(Octant::tsw));
	}
	{
		Octant o = Octant::bse;
		REQUIRE(o.getInteriorNbrOnSide(Side::west) == Octant(Octant::bsw));
		REQUIRE(o.getInteriorNbrOnSide(Side::north) == Octant(Octant::bne));
		REQUIRE(o.getInteriorNbrOnSide(Side::top) == Octant(Octant::tse));
	}
	{
		Octant o = Octant::bnw;
		REQUIRE(o.getInteriorNbrOnSide(Side::east) == Octant(Octant::bne));
		REQUIRE(o.getInteriorNbrOnSide(Side::south) == Octant(Octant::bsw));
		REQUIRE(o.getInteriorNbrOnSide(Side::top) == Octant(Octant::tnw));
	}
	{
		Octant o = Octant::bne;
		REQUIRE(o.getInteriorNbrOnSide(Side::west) == Octant(Octant::bnw));
		REQUIRE(o.getInteriorNbrOnSide(Side::south) == Octant(Octant::bse));
		REQUIRE(o.getInteriorNbrOnSide(Side::top) == Octant(Octant::tne));
	}
	{
		Octant o = Octant::tsw;
		REQUIRE(o.getInteriorNbrOnSide(Side::east) == Octant(Octant::tse));
		REQUIRE(o.getInteriorNbrOnSide(Side::north) == Octant(Octant::tnw));
		REQUIRE(o.getInteriorNbrOnSide(Side::bottom) == Octant(Octant::bsw));
	}
	{
		Octant o = Octant::tse;
		REQUIRE(o.getInteriorNbrOnSide(Side::west) == Octant(Octant::tsw));
		REQUIRE(o.getInteriorNbrOnSide(Side::north) == Octant(Octant::tne));
		REQUIRE(o.getInteriorNbrOnSide(Side::bottom) == Octant(Octant::bse));
	}
	{
		Octant o = Octant::tnw;
		REQUIRE(o.getInteriorNbrOnSide(Side::east) == Octant(Octant::tne));
		REQUIRE(o.getInteriorNbrOnSide(Side::south) == Octant(Octant::tsw));
		REQUIRE(o.getInteriorNbrOnSide(Side::bottom) == Octant(Octant::bnw));
	}
	{
		Octant o = Octant::tne;
		REQUIRE(o.getInteriorNbrOnSide(Side::west) == Octant(Octant::tnw));
		REQUIRE(o.getInteriorNbrOnSide(Side::south) == Octant(Octant::tse));
		REQUIRE(o.getInteriorNbrOnSide(Side::bottom) == Octant(Octant::bne));
	}
}
TEST_CASE("Octant getExteriorNbrOnSide() is as expected", "[Octant]")
{
	{
		Octant o = Octant::bsw;
		REQUIRE(o.getExteriorNbrOnSide(Side::west) == Octant(Octant::bse));
		REQUIRE(o.getExteriorNbrOnSide(Side::south) == Octant(Octant::bnw));
		REQUIRE(o.getExteriorNbrOnSide(Side::bottom) == Octant(Octant::tsw));
	}
	{
		Octant o = Octant::bse;
		REQUIRE(o.getExteriorNbrOnSide(Side::east) == Octant(Octant::bsw));
		REQUIRE(o.getExteriorNbrOnSide(Side::south) == Octant(Octant::bne));
		REQUIRE(o.getExteriorNbrOnSide(Side::bottom) == Octant(Octant::tse));
	}
	{
		Octant o = Octant::bnw;
		REQUIRE(o.getExteriorNbrOnSide(Side::west) == Octant(Octant::bne));
		REQUIRE(o.getExteriorNbrOnSide(Side::north) == Octant(Octant::bsw));
		REQUIRE(o.getExteriorNbrOnSide(Side::bottom) == Octant(Octant::tnw));
	}
	{
		Octant o = Octant::bne;
		REQUIRE(o.getExteriorNbrOnSide(Side::east) == Octant(Octant::bnw));
		REQUIRE(o.getExteriorNbrOnSide(Side::north) == Octant(Octant::bse));
		REQUIRE(o.getExteriorNbrOnSide(Side::bottom) == Octant(Octant::tne));
	}
	{
		Octant o = Octant::tsw;
		REQUIRE(o.getExteriorNbrOnSide(Side::west) == Octant(Octant::tse));
		REQUIRE(o.getExteriorNbrOnSide(Side::south) == Octant(Octant::tnw));
		REQUIRE(o.getExteriorNbrOnSide(Side::top) == Octant(Octant::bsw));
	}
	{
		Octant o = Octant::tse;
		REQUIRE(o.getExteriorNbrOnSide(Side::east) == Octant(Octant::tsw));
		REQUIRE(o.getExteriorNbrOnSide(Side::south) == Octant(Octant::tne));
		REQUIRE(o.getExteriorNbrOnSide(Side::top) == Octant(Octant::bse));
	}
	{
		Octant o = Octant::tnw;
		REQUIRE(o.getExteriorNbrOnSide(Side::west) == Octant(Octant::tne));
		REQUIRE(o.getExteriorNbrOnSide(Side::north) == Octant(Octant::tsw));
		REQUIRE(o.getExteriorNbrOnSide(Side::top) == Octant(Octant::bnw));
	}
	{
		Octant o = Octant::tne;
		REQUIRE(o.getExteriorNbrOnSide(Side::east) == Octant(Octant::tnw));
		REQUIRE(o.getExteriorNbrOnSide(Side::north) == Octant(Octant::tse));
		REQUIRE(o.getExteriorNbrOnSide(Side::top) == Octant(Octant::bne));
	}
}
TEST_CASE("Octant getInteriorSides() is as expected", "[Octant]")
{
	{
		Octant o     = Octant::bsw;
		auto   array = o.getInteriorSides();
		REQUIRE(array[0] == Side(Side::east));
		REQUIRE(array[1] == Side(Side::north));
		REQUIRE(array[2] == Side(Side::top));
	}
	{
		Octant o     = Octant::bse;
		auto   array = o.getInteriorSides();
		REQUIRE(array[0] == Side(Side::west));
		REQUIRE(array[1] == Side(Side::north));
		REQUIRE(array[2] == Side(Side::top));
	}
	{
		Octant o     = Octant::bnw;
		auto   array = o.getInteriorSides();
		REQUIRE(array[0] == Side(Side::east));
		REQUIRE(array[1] == Side(Side::south));
		REQUIRE(array[2] == Side(Side::top));
	}
	{
		Octant o     = Octant::bne;
		auto   array = o.getInteriorSides();
		REQUIRE(array[0] == Side(Side::west));
		REQUIRE(array[1] == Side(Side::south));
		REQUIRE(array[2] == Side(Side::top));
	}
	{
		Octant o     = Octant::tsw;
		auto   array = o.getInteriorSides();
		REQUIRE(array[0] == Side(Side::east));
		REQUIRE(array[1] == Side(Side::north));
		REQUIRE(array[2] == Side(Side::bottom));
	}
	{
		Octant o     = Octant::tse;
		auto   array = o.getInteriorSides();
		REQUIRE(array[0] == Side(Side::west));
		REQUIRE(array[1] == Side(Side::north));
		REQUIRE(array[2] == Side(Side::bottom));
	}
	{
		Octant o     = Octant::tnw;
		auto   array = o.getInteriorSides();
		REQUIRE(array[0] == Side(Side::east));
		REQUIRE(array[1] == Side(Side::south));
		REQUIRE(array[2] == Side(Side::bottom));
	}
	{
		Octant o     = Octant::tne;
		auto   array = o.getInteriorSides();
		REQUIRE(array[0] == Side(Side::west));
		REQUIRE(array[1] == Side(Side::south));
		REQUIRE(array[2] == Side(Side::bottom));
	}
}
TEST_CASE("Octant getExteriorSides() is as expected", "[Octant]")
{
	{
		Octant o     = Octant::bsw;
		auto   array = o.getExteriorSides();
		REQUIRE(array[0] == Side(Side::west));
		REQUIRE(array[1] == Side(Side::south));
		REQUIRE(array[2] == Side(Side::bottom));
	}
	{
		Octant o     = Octant::bse;
		auto   array = o.getExteriorSides();
		REQUIRE(array[0] == Side(Side::east));
		REQUIRE(array[1] == Side(Side::south));
		REQUIRE(array[2] == Side(Side::bottom));
	}
	{
		Octant o     = Octant::bnw;
		auto   array = o.getExteriorSides();
		REQUIRE(array[0] == Side(Side::west));
		REQUIRE(array[1] == Side(Side::north));
		REQUIRE(array[2] == Side(Side::bottom));
	}
	{
		Octant o     = Octant::bne;
		auto   array = o.getExteriorSides();
		REQUIRE(array[0] == Side(Side::east));
		REQUIRE(array[1] == Side(Side::north));
		REQUIRE(array[2] == Side(Side::bottom));
	}
	{
		Octant o     = Octant::tsw;
		auto   array = o.getExteriorSides();
		REQUIRE(array[0] == Side(Side::west));
		REQUIRE(array[1] == Side(Side::south));
		REQUIRE(array[2] == Side(Side::top));
	}
	{
		Octant o     = Octant::tse;
		auto   array = o.getExteriorSides();
		REQUIRE(array[0] == Side(Side::east));
		REQUIRE(array[1] == Side(Side::south));
		REQUIRE(array[2] == Side(Side::top));
	}
	{
		Octant o     = Octant::tnw;
		auto   array = o.getExteriorSides();
		REQUIRE(array[0] == Side(Side::west));
		REQUIRE(array[1] == Side(Side::north));
		REQUIRE(array[2] == Side(Side::top));
	}
	{
		Octant o     = Octant::tne;
		auto   array = o.getExteriorSides();
		REQUIRE(array[0] == Side(Side::east));
		REQUIRE(array[1] == Side(Side::north));
		REQUIRE(array[2] == Side(Side::top));
	}
}
TEST_CASE("Octant isOnSide() is as expected", "[Octant]")
{
	{
		Octant o = Octant::bsw;
		REQUIRE(o.isOnSide(Side::west));
		REQUIRE(!o.isOnSide(Side::east));
		REQUIRE(o.isOnSide(Side::south));
		REQUIRE(!o.isOnSide(Side::north));
		REQUIRE(o.isOnSide(Side::bottom));
		REQUIRE(!o.isOnSide(Side::top));
	}
	{
		Octant o = Octant::bse;
		REQUIRE(!o.isOnSide(Side::west));
		REQUIRE(o.isOnSide(Side::east));
		REQUIRE(o.isOnSide(Side::south));
		REQUIRE(!o.isOnSide(Side::north));
		REQUIRE(o.isOnSide(Side::bottom));
		REQUIRE(!o.isOnSide(Side::top));
	}
	{
		Octant o = Octant::bnw;
		REQUIRE(o.isOnSide(Side::west));
		REQUIRE(!o.isOnSide(Side::east));
		REQUIRE(!o.isOnSide(Side::south));
		REQUIRE(o.isOnSide(Side::north));
		REQUIRE(o.isOnSide(Side::bottom));
		REQUIRE(!o.isOnSide(Side::top));
	}
	{
		Octant o = Octant::bne;
		REQUIRE(!o.isOnSide(Side::west));
		REQUIRE(o.isOnSide(Side::east));
		REQUIRE(!o.isOnSide(Side::south));
		REQUIRE(o.isOnSide(Side::north));
		REQUIRE(o.isOnSide(Side::bottom));
		REQUIRE(!o.isOnSide(Side::top));
	}
	{
		Octant o = Octant::tsw;
		REQUIRE(o.isOnSide(Side::west));
		REQUIRE(!o.isOnSide(Side::east));
		REQUIRE(o.isOnSide(Side::south));
		REQUIRE(!o.isOnSide(Side::north));
		REQUIRE(!o.isOnSide(Side::bottom));
		REQUIRE(o.isOnSide(Side::top));
	}
	{
		Octant o = Octant::tse;
		REQUIRE(!o.isOnSide(Side::west));
		REQUIRE(o.isOnSide(Side::east));
		REQUIRE(o.isOnSide(Side::south));
		REQUIRE(!o.isOnSide(Side::north));
		REQUIRE(!o.isOnSide(Side::bottom));
		REQUIRE(o.isOnSide(Side::top));
	}
	{
		Octant o = Octant::tnw;
		REQUIRE(o.isOnSide(Side::west));
		REQUIRE(!o.isOnSide(Side::east));
		REQUIRE(!o.isOnSide(Side::south));
		REQUIRE(o.isOnSide(Side::north));
		REQUIRE(!o.isOnSide(Side::bottom));
		REQUIRE(o.isOnSide(Side::top));
	}
	{
		Octant o = Octant::tne;
		REQUIRE(!o.isOnSide(Side::west));
		REQUIRE(o.isOnSide(Side::east));
		REQUIRE(!o.isOnSide(Side::south));
		REQUIRE(o.isOnSide(Side::north));
		REQUIRE(!o.isOnSide(Side::bottom));
		REQUIRE(o.isOnSide(Side::top));
	}
}
TEST_CASE("Octant getValues() is as expected", "[Octant]")
{
	std::array<Octant, 8> values = Octant::getValues();
	REQUIRE(values[0] == Octant(Octant::bsw));
	REQUIRE(values[1] == Octant(Octant::bse));
	REQUIRE(values[2] == Octant(Octant::bnw));
	REQUIRE(values[3] == Octant(Octant::bne));
	REQUIRE(values[4] == Octant(Octant::tsw));
	REQUIRE(values[5] == Octant(Octant::tse));
	REQUIRE(values[6] == Octant(Octant::tnw));
	REQUIRE(values[7] == Octant(Octant::tne));
}
TEST_CASE("Octant getValuesOnSide() is as expected", "[Octant]")
{
	{
		std::array<Octant, 4> values = Octant::getValuesOnSide(Side::west);
		REQUIRE(values[0] == Octant(Octant::bsw));
		REQUIRE(values[1] == Octant(Octant::bnw));
		REQUIRE(values[2] == Octant(Octant::tsw));
		REQUIRE(values[3] == Octant(Octant::tnw));
	}
	{
		std::array<Octant, 4> values = Octant::getValuesOnSide(Side::east);
		REQUIRE(values[0] == Octant(Octant::bse));
		REQUIRE(values[1] == Octant(Octant::bne));
		REQUIRE(values[2] == Octant(Octant::tse));
		REQUIRE(values[3] == Octant(Octant::tne));
	}
	{
		std::array<Octant, 4> values = Octant::getValuesOnSide(Side::south);
		REQUIRE(values[0] == Octant(Octant::bsw));
		REQUIRE(values[1] == Octant(Octant::bse));
		REQUIRE(values[2] == Octant(Octant::tsw));
		REQUIRE(values[3] == Octant(Octant::tse));
	}
	{
		std::array<Octant, 4> values = Octant::getValuesOnSide(Side::north);
		REQUIRE(values[0] == Octant(Octant::bnw));
		REQUIRE(values[1] == Octant(Octant::bne));
		REQUIRE(values[2] == Octant(Octant::tnw));
		REQUIRE(values[3] == Octant(Octant::tne));
	}
	{
		std::array<Octant, 4> values = Octant::getValuesOnSide(Side::bottom);
		REQUIRE(values[0] == Octant(Octant::bsw));
		REQUIRE(values[1] == Octant(Octant::bse));
		REQUIRE(values[2] == Octant(Octant::bnw));
		REQUIRE(values[3] == Octant(Octant::bne));
	}
	{
		std::array<Octant, 4> values = Octant::getValuesOnSide(Side::top);
		REQUIRE(values[0] == Octant(Octant::tsw));
		REQUIRE(values[1] == Octant(Octant::tse));
		REQUIRE(values[2] == Octant(Octant::tnw));
		REQUIRE(values[3] == Octant(Octant::tne));
	}
}
