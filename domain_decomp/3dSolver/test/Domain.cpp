#include "../Domain.h"
#include "catch.hpp"
using namespace std;
TEST_CASE("NormalNbrInfo getNbrType works", "[Domain]")
{
	NbrInfo *info = new NormalNbrInfo();
	REQUIRE(info->getNbrType() == NbrType::Normal);
	delete info;
}

TEST_CASE("NormalNbrInfo Serialization/Deserialization", "[Domain]")
{
	NormalNbrInfo info;
	info.id   = 5;
	info.rank = 1;
	// serialize and then deserialize
	char *buff = new char[info.serialize(nullptr)];
	info.serialize(buff);
	NormalNbrInfo out;
	out.deserialize(buff);
	delete[] buff;
	REQUIRE(out.id == 5);
	REQUIRE(out.rank == 1);
}
TEST_CASE("CoarseNbrInfo Serialization/Deserialization", "[Domain]")
{
	CoarseNbrInfo info;
	info.id             = 5;
	info.rank           = 1;
	info.quad_on_coarse = 2;
	// serialize and then deserialize
	char *buff = new char[info.serialize(nullptr)];
	info.serialize(buff);
	CoarseNbrInfo out;
	out.deserialize(buff);
	delete[] buff;
	REQUIRE(out.id == 5);
	REQUIRE(out.rank == 1);
	REQUIRE(out.quad_on_coarse == 2);
}
TEST_CASE("FineNbrInfo Serialization/Deserialization", "[Domain]")
{
	FineNbrInfo info;
	info.ids[0]   = 1;
	info.ids[1]   = 2;
	info.ids[2]   = 3;
	info.ids[3]   = 4;
	info.ranks[0] = 9;
	info.ranks[1] = 8;
	info.ranks[2] = 7;
	info.ranks[3] = 6;
	// serialize and then deserialize
	char *buff = new char[info.serialize(nullptr)];
	info.serialize(buff);
	FineNbrInfo out;
	out.deserialize(buff);
	delete[] buff;
	REQUIRE(out.ids[0] == 1);
	REQUIRE(out.ids[1] == 2);
	REQUIRE(out.ids[2] == 3);
	REQUIRE(out.ids[3] == 4);
	REQUIRE(out.ranks[0] == 9);
	REQUIRE(out.ranks[1] == 8);
	REQUIRE(out.ranks[2] == 7);
	REQUIRE(out.ranks[3] == 6);
}
TEST_CASE("Domain Serialization/Deserialization", "[Domain]")
{
	Domain *d_ptr                = new Domain;
	Domain &d                    = *d_ptr;
	d.id                         = 0;
	d.getNbrInfoPtr(Side::north) = new NormalNbrInfo(1);
	d.getNbrInfoPtr(Side::east)  = new CoarseNbrInfo(2, 3);
	d.getNbrInfoPtr(Side::south) = new FineNbrInfo({3, 4, 5, 6});

	// serialize and then deserialize
	char *buff = new char[d.serialize(nullptr)];
	d.serialize(buff);
	delete d_ptr;
	Domain out;
	out.deserialize(buff);
	delete[] buff;

	// check that deserialized version has the same information
	REQUIRE(out.id == 0);

	REQUIRE(!out.hasNbr(Side::west));

	REQUIRE(out.hasNbr(Side::east));
	REQUIRE(out.getNbrType(Side::east) == NbrType::Coarse);
	REQUIRE(out.getCoarseNbrInfo(Side::east).id == 2);
	REQUIRE(out.getCoarseNbrInfo(Side::east).quad_on_coarse == 3);

	REQUIRE(out.hasNbr(Side::south));
	REQUIRE(out.getNbrType(Side::south) == NbrType::Fine);
	REQUIRE(out.getFineNbrInfo(Side::south).ids[0] == 3);
	REQUIRE(out.getFineNbrInfo(Side::south).ids[1] == 4);
	REQUIRE(out.getFineNbrInfo(Side::south).ids[2] == 5);
	REQUIRE(out.getFineNbrInfo(Side::south).ids[3] == 6);

	REQUIRE(out.hasNbr(Side::north));
	REQUIRE(out.getNbrType(Side::north) == NbrType::Normal);
	REQUIRE(out.getNormalNbrInfo(Side::north).id == 1);

	REQUIRE(!out.hasNbr(Side::bottom));
	REQUIRE(!out.hasNbr(Side::top));
}
