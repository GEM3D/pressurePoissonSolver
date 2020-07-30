#include <Thunderegg/Side.h>
#include "catch.hpp"
using namespace std;
TEST_CASE("Orthant<3> Default constructor works", "[Octant]")
{
	Orthant<3> o;
	CHECK(o.toInt() == -1);
}
TEST_CASE("Orthant<3> int constructor and toInt() works", "[Octant]")
{
	{
		Orthant<3> o(5);
		CHECK(o.toInt() == 5);
	}
	{
		Orthant<3> o(-8);
		CHECK(o.toInt() == -8);
	}
}
TEST_CASE("Otant == operator works", "[Octant]")
{
	{
		Orthant<3> o(5);
		CHECK(!(o == 4));
		CHECK(o == 5);
		CHECK(!(o == 6));
	}
	{
		Side<3> o(0);
		CHECK(!(o == -1));
		CHECK(o == 0);
		CHECK(!(o == 1));
	}
}
TEST_CASE("Orthant<3> values are as expected", "[Octant]")
{
	{
		Orthant<3> o(0);
		Orthant<3> expected = Orthant<3>::bsw;
		CHECK(o == expected);
	}
	{
		Orthant<3> o(1);
		Orthant<3> expected = Orthant<3>::bse;
		CHECK(o == expected);
	}
	{
		Orthant<3> o(2);
		Orthant<3> expected = Orthant<3>::bnw;
		CHECK(o == expected);
	}
	{
		Orthant<3> o(3);
		Orthant<3> expected = Orthant<3>::bne;
		CHECK(o == expected);
	}
	{
		Orthant<3> o(4);
		Orthant<3> expected = Orthant<3>::tsw;
		CHECK(o == expected);
	}
	{
		Orthant<3> o(5);
		Orthant<3> expected = Orthant<3>::tse;
		CHECK(o == expected);
	}
	{
		Orthant<3> o(6);
		Orthant<3> expected = Orthant<3>::tnw;
		CHECK(o == expected);
	}
	{
		Orthant<3> o(7);
		Orthant<3> expected = Orthant<3>::tne;
		CHECK(o == expected);
	}
}
TEST_CASE("Orthant<3> getInteriorNbrOnSide() is as expected", "[Octant]")
{
	{
		Orthant<3> o = Orthant<3>::bsw;
		CHECK(o.getInteriorNbrOnSide(Side<3>::east) == Orthant<3>(Orthant<3>::bse));
		CHECK(o.getInteriorNbrOnSide(Side<3>::north) == Orthant<3>(Orthant<3>::bnw));
		CHECK(o.getInteriorNbrOnSide(Side<3>::top) == Orthant<3>(Orthant<3>::tsw));
	}
	{
		Orthant<3> o = Orthant<3>::bse;
		CHECK(o.getInteriorNbrOnSide(Side<3>::west) == Orthant<3>(Orthant<3>::bsw));
		CHECK(o.getInteriorNbrOnSide(Side<3>::north) == Orthant<3>(Orthant<3>::bne));
		CHECK(o.getInteriorNbrOnSide(Side<3>::top) == Orthant<3>(Orthant<3>::tse));
	}
	{
		Orthant<3> o = Orthant<3>::bnw;
		CHECK(o.getInteriorNbrOnSide(Side<3>::east) == Orthant<3>(Orthant<3>::bne));
		CHECK(o.getInteriorNbrOnSide(Side<3>::south) == Orthant<3>(Orthant<3>::bsw));
		CHECK(o.getInteriorNbrOnSide(Side<3>::top) == Orthant<3>(Orthant<3>::tnw));
	}
	{
		Orthant<3> o = Orthant<3>::bne;
		CHECK(o.getInteriorNbrOnSide(Side<3>::west) == Orthant<3>(Orthant<3>::bnw));
		CHECK(o.getInteriorNbrOnSide(Side<3>::south) == Orthant<3>(Orthant<3>::bse));
		CHECK(o.getInteriorNbrOnSide(Side<3>::top) == Orthant<3>(Orthant<3>::tne));
	}
	{
		Orthant<3> o = Orthant<3>::tsw;
		CHECK(o.getInteriorNbrOnSide(Side<3>::east) == Orthant<3>(Orthant<3>::tse));
		CHECK(o.getInteriorNbrOnSide(Side<3>::north) == Orthant<3>(Orthant<3>::tnw));
		CHECK(o.getInteriorNbrOnSide(Side<3>::bottom) == Orthant<3>(Orthant<3>::bsw));
	}
	{
		Orthant<3> o = Orthant<3>::tse;
		CHECK(o.getInteriorNbrOnSide(Side<3>::west) == Orthant<3>(Orthant<3>::tsw));
		CHECK(o.getInteriorNbrOnSide(Side<3>::north) == Orthant<3>(Orthant<3>::tne));
		CHECK(o.getInteriorNbrOnSide(Side<3>::bottom) == Orthant<3>(Orthant<3>::bse));
	}
	{
		Orthant<3> o = Orthant<3>::tnw;
		CHECK(o.getInteriorNbrOnSide(Side<3>::east) == Orthant<3>(Orthant<3>::tne));
		CHECK(o.getInteriorNbrOnSide(Side<3>::south) == Orthant<3>(Orthant<3>::tsw));
		CHECK(o.getInteriorNbrOnSide(Side<3>::bottom) == Orthant<3>(Orthant<3>::bnw));
	}
	{
		Orthant<3> o = Orthant<3>::tne;
		CHECK(o.getInteriorNbrOnSide(Side<3>::west) == Orthant<3>(Orthant<3>::tnw));
		CHECK(o.getInteriorNbrOnSide(Side<3>::south) == Orthant<3>(Orthant<3>::tse));
		CHECK(o.getInteriorNbrOnSide(Side<3>::bottom) == Orthant<3>(Orthant<3>::bne));
	}
}
TEST_CASE("Orthant<3> getExteriorNbrOnSide() is as expected", "[Octant]")
{
	{
		Orthant<3> o = Orthant<3>::bsw;
		CHECK(o.getExteriorNbrOnSide(Side<3>::west) == Orthant<3>(Orthant<3>::bse));
		CHECK(o.getExteriorNbrOnSide(Side<3>::south) == Orthant<3>(Orthant<3>::bnw));
		CHECK(o.getExteriorNbrOnSide(Side<3>::bottom) == Orthant<3>(Orthant<3>::tsw));
	}
	{
		Orthant<3> o = Orthant<3>::bse;
		CHECK(o.getExteriorNbrOnSide(Side<3>::east) == Orthant<3>(Orthant<3>::bsw));
		CHECK(o.getExteriorNbrOnSide(Side<3>::south) == Orthant<3>(Orthant<3>::bne));
		CHECK(o.getExteriorNbrOnSide(Side<3>::bottom) == Orthant<3>(Orthant<3>::tse));
	}
	{
		Orthant<3> o = Orthant<3>::bnw;
		CHECK(o.getExteriorNbrOnSide(Side<3>::west) == Orthant<3>(Orthant<3>::bne));
		CHECK(o.getExteriorNbrOnSide(Side<3>::north) == Orthant<3>(Orthant<3>::bsw));
		CHECK(o.getExteriorNbrOnSide(Side<3>::bottom) == Orthant<3>(Orthant<3>::tnw));
	}
	{
		Orthant<3> o = Orthant<3>::bne;
		CHECK(o.getExteriorNbrOnSide(Side<3>::east) == Orthant<3>(Orthant<3>::bnw));
		CHECK(o.getExteriorNbrOnSide(Side<3>::north) == Orthant<3>(Orthant<3>::bse));
		CHECK(o.getExteriorNbrOnSide(Side<3>::bottom) == Orthant<3>(Orthant<3>::tne));
	}
	{
		Orthant<3> o = Orthant<3>::tsw;
		CHECK(o.getExteriorNbrOnSide(Side<3>::west) == Orthant<3>(Orthant<3>::tse));
		CHECK(o.getExteriorNbrOnSide(Side<3>::south) == Orthant<3>(Orthant<3>::tnw));
		CHECK(o.getExteriorNbrOnSide(Side<3>::top) == Orthant<3>(Orthant<3>::bsw));
	}
	{
		Orthant<3> o = Orthant<3>::tse;
		CHECK(o.getExteriorNbrOnSide(Side<3>::east) == Orthant<3>(Orthant<3>::tsw));
		CHECK(o.getExteriorNbrOnSide(Side<3>::south) == Orthant<3>(Orthant<3>::tne));
		CHECK(o.getExteriorNbrOnSide(Side<3>::top) == Orthant<3>(Orthant<3>::bse));
	}
	{
		Orthant<3> o = Orthant<3>::tnw;
		CHECK(o.getExteriorNbrOnSide(Side<3>::west) == Orthant<3>(Orthant<3>::tne));
		CHECK(o.getExteriorNbrOnSide(Side<3>::north) == Orthant<3>(Orthant<3>::tsw));
		CHECK(o.getExteriorNbrOnSide(Side<3>::top) == Orthant<3>(Orthant<3>::bnw));
	}
	{
		Orthant<3> o = Orthant<3>::tne;
		CHECK(o.getExteriorNbrOnSide(Side<3>::east) == Orthant<3>(Orthant<3>::tnw));
		CHECK(o.getExteriorNbrOnSide(Side<3>::north) == Orthant<3>(Orthant<3>::tse));
		CHECK(o.getExteriorNbrOnSide(Side<3>::top) == Orthant<3>(Orthant<3>::bne));
	}
}
TEST_CASE("Orthant<3> getInteriorSides() is as expected", "[Octant]")
{
	{
		Orthant<3> o     = Orthant<3>::bsw;
		auto       array = o.getInteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::east));
		CHECK(array[1] == Side<3>(Side<3>::north));
		CHECK(array[2] == Side<3>(Side<3>::top));
	}
	{
		Orthant<3> o     = Orthant<3>::bse;
		auto       array = o.getInteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::west));
		CHECK(array[1] == Side<3>(Side<3>::north));
		CHECK(array[2] == Side<3>(Side<3>::top));
	}
	{
		Orthant<3> o     = Orthant<3>::bnw;
		auto       array = o.getInteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::east));
		CHECK(array[1] == Side<3>(Side<3>::south));
		CHECK(array[2] == Side<3>(Side<3>::top));
	}
	{
		Orthant<3> o     = Orthant<3>::bne;
		auto       array = o.getInteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::west));
		CHECK(array[1] == Side<3>(Side<3>::south));
		CHECK(array[2] == Side<3>(Side<3>::top));
	}
	{
		Orthant<3> o     = Orthant<3>::tsw;
		auto       array = o.getInteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::east));
		CHECK(array[1] == Side<3>(Side<3>::north));
		CHECK(array[2] == Side<3>(Side<3>::bottom));
	}
	{
		Orthant<3> o     = Orthant<3>::tse;
		auto       array = o.getInteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::west));
		CHECK(array[1] == Side<3>(Side<3>::north));
		CHECK(array[2] == Side<3>(Side<3>::bottom));
	}
	{
		Orthant<3> o     = Orthant<3>::tnw;
		auto       array = o.getInteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::east));
		CHECK(array[1] == Side<3>(Side<3>::south));
		CHECK(array[2] == Side<3>(Side<3>::bottom));
	}
	{
		Orthant<3> o     = Orthant<3>::tne;
		auto       array = o.getInteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::west));
		CHECK(array[1] == Side<3>(Side<3>::south));
		CHECK(array[2] == Side<3>(Side<3>::bottom));
	}
}
TEST_CASE("Orthant<3> getExteriorSides() is as expected", "[Octant]")
{
	{
		Orthant<3> o     = Orthant<3>::bsw;
		auto       array = o.getExteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::west));
		CHECK(array[1] == Side<3>(Side<3>::south));
		CHECK(array[2] == Side<3>(Side<3>::bottom));
	}
	{
		Orthant<3> o     = Orthant<3>::bse;
		auto       array = o.getExteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::east));
		CHECK(array[1] == Side<3>(Side<3>::south));
		CHECK(array[2] == Side<3>(Side<3>::bottom));
	}
	{
		Orthant<3> o     = Orthant<3>::bnw;
		auto       array = o.getExteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::west));
		CHECK(array[1] == Side<3>(Side<3>::north));
		CHECK(array[2] == Side<3>(Side<3>::bottom));
	}
	{
		Orthant<3> o     = Orthant<3>::bne;
		auto       array = o.getExteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::east));
		CHECK(array[1] == Side<3>(Side<3>::north));
		CHECK(array[2] == Side<3>(Side<3>::bottom));
	}
	{
		Orthant<3> o     = Orthant<3>::tsw;
		auto       array = o.getExteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::west));
		CHECK(array[1] == Side<3>(Side<3>::south));
		CHECK(array[2] == Side<3>(Side<3>::top));
	}
	{
		Orthant<3> o     = Orthant<3>::tse;
		auto       array = o.getExteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::east));
		CHECK(array[1] == Side<3>(Side<3>::south));
		CHECK(array[2] == Side<3>(Side<3>::top));
	}
	{
		Orthant<3> o     = Orthant<3>::tnw;
		auto       array = o.getExteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::west));
		CHECK(array[1] == Side<3>(Side<3>::north));
		CHECK(array[2] == Side<3>(Side<3>::top));
	}
	{
		Orthant<3> o     = Orthant<3>::tne;
		auto       array = o.getExteriorSides();
		CHECK(array[0] == Side<3>(Side<3>::east));
		CHECK(array[1] == Side<3>(Side<3>::north));
		CHECK(array[2] == Side<3>(Side<3>::top));
	}
}
TEST_CASE("Orthant<3> isOnSide() is as expected", "[Octant]")
{
	{
		Orthant<3> o = Orthant<3>::bsw;
		CHECK(o.isOnSide(Side<3>::west));
		CHECK(!o.isOnSide(Side<3>::east));
		CHECK(o.isOnSide(Side<3>::south));
		CHECK(!o.isOnSide(Side<3>::north));
		CHECK(o.isOnSide(Side<3>::bottom));
		CHECK(!o.isOnSide(Side<3>::top));
	}
	{
		Orthant<3> o = Orthant<3>::bse;
		CHECK(!o.isOnSide(Side<3>::west));
		CHECK(o.isOnSide(Side<3>::east));
		CHECK(o.isOnSide(Side<3>::south));
		CHECK(!o.isOnSide(Side<3>::north));
		CHECK(o.isOnSide(Side<3>::bottom));
		CHECK(!o.isOnSide(Side<3>::top));
	}
	{
		Orthant<3> o = Orthant<3>::bnw;
		CHECK(o.isOnSide(Side<3>::west));
		CHECK(!o.isOnSide(Side<3>::east));
		CHECK(!o.isOnSide(Side<3>::south));
		CHECK(o.isOnSide(Side<3>::north));
		CHECK(o.isOnSide(Side<3>::bottom));
		CHECK(!o.isOnSide(Side<3>::top));
	}
	{
		Orthant<3> o = Orthant<3>::bne;
		CHECK(!o.isOnSide(Side<3>::west));
		CHECK(o.isOnSide(Side<3>::east));
		CHECK(!o.isOnSide(Side<3>::south));
		CHECK(o.isOnSide(Side<3>::north));
		CHECK(o.isOnSide(Side<3>::bottom));
		CHECK(!o.isOnSide(Side<3>::top));
	}
	{
		Orthant<3> o = Orthant<3>::tsw;
		CHECK(o.isOnSide(Side<3>::west));
		CHECK(!o.isOnSide(Side<3>::east));
		CHECK(o.isOnSide(Side<3>::south));
		CHECK(!o.isOnSide(Side<3>::north));
		CHECK(!o.isOnSide(Side<3>::bottom));
		CHECK(o.isOnSide(Side<3>::top));
	}
	{
		Orthant<3> o = Orthant<3>::tse;
		CHECK(!o.isOnSide(Side<3>::west));
		CHECK(o.isOnSide(Side<3>::east));
		CHECK(o.isOnSide(Side<3>::south));
		CHECK(!o.isOnSide(Side<3>::north));
		CHECK(!o.isOnSide(Side<3>::bottom));
		CHECK(o.isOnSide(Side<3>::top));
	}
	{
		Orthant<3> o = Orthant<3>::tnw;
		CHECK(o.isOnSide(Side<3>::west));
		CHECK(!o.isOnSide(Side<3>::east));
		CHECK(!o.isOnSide(Side<3>::south));
		CHECK(o.isOnSide(Side<3>::north));
		CHECK(!o.isOnSide(Side<3>::bottom));
		CHECK(o.isOnSide(Side<3>::top));
	}
	{
		Orthant<3> o = Orthant<3>::tne;
		CHECK(!o.isOnSide(Side<3>::west));
		CHECK(o.isOnSide(Side<3>::east));
		CHECK(!o.isOnSide(Side<3>::south));
		CHECK(o.isOnSide(Side<3>::north));
		CHECK(!o.isOnSide(Side<3>::bottom));
		CHECK(o.isOnSide(Side<3>::top));
	}
}
TEST_CASE("Orthant<3> getValues() is as expected", "[Octant]")
{
	std::array<Orthant<3>, 8> values = Orthant<3>::getValues();
	CHECK(values[0] == Orthant<3>(Orthant<3>::bsw));
	CHECK(values[1] == Orthant<3>(Orthant<3>::bse));
	CHECK(values[2] == Orthant<3>(Orthant<3>::bnw));
	CHECK(values[3] == Orthant<3>(Orthant<3>::bne));
	CHECK(values[4] == Orthant<3>(Orthant<3>::tsw));
	CHECK(values[5] == Orthant<3>(Orthant<3>::tse));
	CHECK(values[6] == Orthant<3>(Orthant<3>::tnw));
	CHECK(values[7] == Orthant<3>(Orthant<3>::tne));
}
TEST_CASE("Orthant<3> getValuesOnSide() is as expected", "[Octant]")
{
	SECTION("Side<3>::west")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::west);
		CHECK(values[0] == Orthant<3>(Orthant<3>::bsw));
		CHECK(values[1] == Orthant<3>(Orthant<3>::bnw));
		CHECK(values[2] == Orthant<3>(Orthant<3>::tsw));
		CHECK(values[3] == Orthant<3>(Orthant<3>::tnw));
	}
	SECTION("Side<3>::east")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::east);
		CHECK(values[0] == Orthant<3>(Orthant<3>::bse));
		CHECK(values[1] == Orthant<3>(Orthant<3>::bne));
		CHECK(values[2] == Orthant<3>(Orthant<3>::tse));
		CHECK(values[3] == Orthant<3>(Orthant<3>::tne));
	}
	SECTION("Side<3>::south")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::south);
		CHECK(values[0] == Orthant<3>(Orthant<3>::bsw));
		CHECK(values[1] == Orthant<3>(Orthant<3>::bse));
		CHECK(values[2] == Orthant<3>(Orthant<3>::tsw));
		CHECK(values[3] == Orthant<3>(Orthant<3>::tse));
	}
	SECTION("Side<3>::north")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::north);
		CHECK(values[0] == Orthant<3>(Orthant<3>::bnw));
		CHECK(values[1] == Orthant<3>(Orthant<3>::bne));
		CHECK(values[2] == Orthant<3>(Orthant<3>::tnw));
		CHECK(values[3] == Orthant<3>(Orthant<3>::tne));
	}
	SECTION("Side<3>::bottom")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::bottom);
		CHECK(values[0] == Orthant<3>(Orthant<3>::bsw));
		CHECK(values[1] == Orthant<3>(Orthant<3>::bse));
		CHECK(values[2] == Orthant<3>(Orthant<3>::bnw));
		CHECK(values[3] == Orthant<3>(Orthant<3>::bne));
	}
	SECTION("Side<3>::top")
	{
		std::array<Orthant<3>, 4> values = Orthant<3>::getValuesOnSide(Side<3>::top);
		CHECK(values[0] == Orthant<3>(Orthant<3>::tsw));
		CHECK(values[1] == Orthant<3>(Orthant<3>::tse));
		CHECK(values[2] == Orthant<3>(Orthant<3>::tnw));
		CHECK(values[3] == Orthant<3>(Orthant<3>::tne));
	}
}
