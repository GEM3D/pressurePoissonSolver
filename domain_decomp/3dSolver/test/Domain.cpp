#include "../Domain.h"
#include "catch.hpp"
TEST_CASE("NormalNbrInfo getNbrType works", "[Domain]")
{
	NbrInfo *info = new NormalNbrInfo();
	REQUIRE(info->getNbrType() == NbrType::Normal);
    delete info;
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
	Domain out = Domain::deserialize(buff);
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
