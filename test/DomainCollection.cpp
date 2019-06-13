#include "../Domain.h"
#include "catch.hpp"
using namespace std;
#if 0
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
			Domain<3> &d = *p.second;
			for (Side<3> s : Side<3>::getValues()) {
				if (d.hasNbr(s)) {
					REQUIRE(d.getNbrType(s) == NbrType::Normal);
					NormalNbrInfo<3> info = d.getNormalNbrInfo(s);
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
			Domain<3> &d = *p.second;
			for (Side<3> s : Side<3>::getValues()) {
				if (d.hasNbr(s)) {
					REQUIRE(d.getNbrType(s) == NbrType::Normal);
					NormalNbrInfo<3> info = d.getNormalNbrInfo(s);
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
			Domain<3> &d = *p.second;
			for (Side<3> s : Side<3>::getValues()) {
				if (d.hasNbr(s)) {
					switch (d.getNbrType(s)) {
						case NbrType::Normal: {
							NormalNbrInfo<3> info = d.getNormalNbrInfo(s);
							REQUIRE(domains[info.id]->getNormalNbrInfo(s.opposite()).id == d.id);
						} break;
						case NbrType::Fine: {
							FineNbrInfo<3> info = d.getFineNbrInfo(s);
							for (int i = 0; i < 4; i++) {
								REQUIRE(domains[info.ids[i]]->getCoarseNbrInfo(s.opposite()).id
								        == d.id);
							}
						} break;
						case NbrType::Coarse: {
							CoarseNbrInfo<3> info = d.getCoarseNbrInfo(s);
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
		for (auto &p : domains) {
			Domain<3> &d = *p.second;
			if (d.refine_level < 3) { REQUIRE(d.id == d.parent_id); }
			// check that neighbor info makes sense
			for (Side<3> s : Side<3>::getValues()) {
				if (d.hasNbr(s)) {
					switch (d.getNbrType(s)) {
						case NbrType::Normal: {
							NormalNbrInfo<3> info = d.getNormalNbrInfo(s);
							REQUIRE(domains[info.id]->getNormalNbrInfo(s.opposite()).id == d.id);
						} break;
						case NbrType::Fine: {
							FineNbrInfo<3> info = d.getFineNbrInfo(s);
							for (int i = 0; i < 4; i++) {
								REQUIRE(domains[info.ids[i]]->getCoarseNbrInfo(s.opposite()).id
								        == d.id);
							}
						} break;
						case NbrType::Coarse: {
							CoarseNbrInfo<3> info = d.getCoarseNbrInfo(s);
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
#endif
