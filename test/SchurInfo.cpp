#include <Thunderegg/SchurInfo.h>
#include "catch.hpp"

using namespace std;
TEST_CASE("SchurInfo enumerateIfaces function works", "[SchurInfo]")
{
	// create schur domain with neighbors only on one proc
	shared_ptr<PatchInfo<3>> d(new PatchInfo<3>);
	d->id                            = 0;
	d->nbr_info[Side<3>::north].reset( new NormalNbrInfo<3>(1));
	d->nbr_info[Side<3>::east].reset(new CoarseNbrInfo<3>(2, 3));
	d->nbr_info[Side<3>::south].reset(new FineNbrInfo<3>({4, 5, 6, 7}));
	SchurInfo<3> sinfo(d);

	map<int,IfaceSet<3>> ifaces;
	map<int,map<int,IfaceSet<3>>> opi;
	set<int> in_procs;
	sinfo.enumerateIfaces(ifaces,opi,in_procs);
	REQUIRE(in_procs.size()==0);
	REQUIRE(opi.size()==0);

}
