#include "../DomainCollection.h"
#include "catch.hpp"
using namespace std;
TEST_CASE("DomainCollection constructors work", "[DomainCollection]")
{
	PetscInitialize(nullptr, nullptr, nullptr, nullptr);
	OctTree uniform("3uni.bin");
	OctTree refined("2refine.bin");
	SECTION("Leaf constructor works on uniform mesh")
	{
		// generate simple tree
		DomainCollection dc(uniform, 10);
		auto             domains = dc.getDomainMap();
		// check that neighbor info makes sense
		for (auto &p : domains) {
			Domain &d = *p.second;
			for (Side s : Side::getValues()) {
				if (d.hasNbr(s)) {
					REQUIRE(d.getNbrType(s) == NbrType::Normal);
					NormalNbrInfo info = d.getNormalNbrInfo(s);
					REQUIRE(domains[info.id]->getNormalNbrInfo(s.opposite()).id == d.id);
				}
			}
		}
	}
	SECTION("Level constructor works on uniform mesh")
	{
		// generate simple tree
		DomainCollection dc(uniform, 2, 10);
		auto             domains = dc.getDomainMap();
		// check that neighbor info makes sense
		for (auto &p : domains) {
			Domain &d = *p.second;
			for (Side s : Side::getValues()) {
				if (d.hasNbr(s)) {
					REQUIRE(d.getNbrType(s) == NbrType::Normal);
					NormalNbrInfo info = d.getNormalNbrInfo(s);
					REQUIRE(domains[info.id]->getNormalNbrInfo(s.opposite()).id == d.id);
				}
			}
		}
	}
	SECTION("Leaf constructor works on refined mesh")
	{
		// generate simple tree
		DomainCollection dc(refined, 10);
		auto             domains = dc.getDomainMap();
		// check that neighbor info makes sense
		for (auto &p : domains) {
			Domain &d = *p.second;
			for (Side s : Side::getValues()) {
				if (d.hasNbr(s)) {
					switch (d.getNbrType(s)) {
						case NbrType::Normal: {
							NormalNbrInfo info = d.getNormalNbrInfo(s);
							REQUIRE(domains[info.id]->getNormalNbrInfo(s.opposite()).id == d.id);
						} break;
						case NbrType::Fine: {
							FineNbrInfo info = d.getFineNbrInfo(s);
							for (int i = 0; i < 4; i++) {
								REQUIRE(domains[info.ids[i]]->getCoarseNbrInfo(s.opposite()).id
								        == d.id);
							}
						} break;
						case NbrType::Coarse: {
							CoarseNbrInfo info = d.getCoarseNbrInfo(s);
							REQUIRE(
							domains[info.id]->getFineNbrInfo(s.opposite()).ids[info.quad_on_coarse]
							== d.id);
						} break;
					}
				}
			}
		}
	}
	SECTION("Level constructor works on refined mesh")
	{
		// generate simple tree
		DomainCollection dc(refined, 3, 10);
		auto             domains = dc.getDomainMap();
		// check that neighbor info makes sense
		for (auto &p : domains) {
			Domain &d = *p.second;
			for (Side s : Side::getValues()) {
				if (d.hasNbr(s)) {
					switch (d.getNbrType(s)) {
						case NbrType::Normal: {
							NormalNbrInfo info = d.getNormalNbrInfo(s);
							REQUIRE(domains[info.id]->getNormalNbrInfo(s.opposite()).id == d.id);
						} break;
						case NbrType::Fine: {
							FineNbrInfo info = d.getFineNbrInfo(s);
							for (int i = 0; i < 4; i++) {
								REQUIRE(domains[info.ids[i]]->getCoarseNbrInfo(s.opposite()).id
								        == d.id);
							}
						} break;
						case NbrType::Coarse: {
							CoarseNbrInfo info = d.getCoarseNbrInfo(s);
							REQUIRE(
							domains[info.id]->getFineNbrInfo(s.opposite()).ids[info.quad_on_coarse]
							== d.id);
						} break;
					}
				}
			}
		}
	}
}
