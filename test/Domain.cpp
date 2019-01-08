#include "../Domain.h"
#include "catch.hpp"
using namespace std;
TEST_CASE("NormalNbrInfo getNbrType works", "[Domain]")
{
	NbrInfo<3> *info = new NormalNbrInfo<3>();
	REQUIRE(info->getNbrType() == NbrType::Normal);
	delete info;
}

TEST_CASE("NormalNbrInfo Serialization/Deserialization", "[Domain]")
{
	NormalNbrInfo<3> info;
	info.id   = 5;
	info.rank = 1;
	// serialize and then deserialize
	char *buff = new char[info.serialize(nullptr)];
	info.serialize(buff);
	NormalNbrInfo<3> out;
	out.deserialize(buff);
	delete[] buff;
	REQUIRE(out.id == 5);
	REQUIRE(out.rank == 1);
}
TEST_CASE("CoarseNbrInfo Serialization/Deserialization", "[Domain]")
{
	CoarseNbrInfo<3> info;
	info.id             = 5;
	info.rank           = 1;
	info.quad_on_coarse = 2;
	// serialize and then deserialize
	char *buff = new char[info.serialize(nullptr)];
	info.serialize(buff);
	CoarseNbrInfo<3> out;
	out.deserialize(buff);
	delete[] buff;
	REQUIRE(out.id == 5);
	REQUIRE(out.rank == 1);
	REQUIRE(out.quad_on_coarse == 2);
}
TEST_CASE("FineNbrInfo Serialization/Deserialization", "[Domain]")
{
	FineNbrInfo<3> info;
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
	FineNbrInfo<3> out;
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
	Domain<3> *d_ptr                = new Domain<3>;
	Domain<3> &d                    = *d_ptr;
	d.id                            = 0;
	d.getNbrInfoPtr(Side<3>::north) = new NormalNbrInfo<3>(1);
	d.getNbrInfoPtr(Side<3>::east)  = new CoarseNbrInfo<3>(2, 3);
	d.getNbrInfoPtr(Side<3>::south) = new FineNbrInfo<3>({3, 4, 5, 6});

	// serialize and then deserialize
	char *buff = new char[d.serialize(nullptr)];
	d.serialize(buff);
	delete d_ptr;
	Domain<3> out;
	out.deserialize(buff);
	delete[] buff;

	// check that deserialized version has the same information
	REQUIRE(out.id == 0);

	REQUIRE(!out.hasNbr(Side<3>::west));

	REQUIRE(out.hasNbr(Side<3>::east));
	REQUIRE(out.getNbrType(Side<3>::east) == NbrType::Coarse);
	REQUIRE(out.getCoarseNbrInfo(Side<3>::east).id == 2);
	REQUIRE(out.getCoarseNbrInfo(Side<3>::east).quad_on_coarse == 3);

	REQUIRE(out.hasNbr(Side<3>::south));
	REQUIRE(out.getNbrType(Side<3>::south) == NbrType::Fine);
	REQUIRE(out.getFineNbrInfo(Side<3>::south).ids[0] == 3);
	REQUIRE(out.getFineNbrInfo(Side<3>::south).ids[1] == 4);
	REQUIRE(out.getFineNbrInfo(Side<3>::south).ids[2] == 5);
	REQUIRE(out.getFineNbrInfo(Side<3>::south).ids[3] == 6);

	REQUIRE(out.hasNbr(Side<3>::north));
	REQUIRE(out.getNbrType(Side<3>::north) == NbrType::Normal);
	REQUIRE(out.getNormalNbrInfo(Side<3>::north).id == 1);

	REQUIRE(!out.hasNbr(Side<3>::bottom));
	REQUIRE(!out.hasNbr(Side<3>::top));
}
