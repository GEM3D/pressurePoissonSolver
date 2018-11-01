#include "Helper.h"
#include "AvgRstr.h"
#include "DrctIntp.h"
#include "FFTBlockJacobiSmoother.h"
#include "MatOp.h"
#include "MatrixHelper.h"
#include "TriLinIntp.h"
#include "VCycle.h"
#include "WCycle.h"
#include "WrapOp.h"
#include <fstream>
using namespace std;
using namespace GMG;
using nlohmann::json;
Helper::Helper(int n, OctTree t, std::vector<std::shared_ptr<DomainCollection<3>>> dcs,
               std::shared_ptr<SchurHelper<3>> sh, std::string config_file)
{
	ifstream config_stream(config_file);
	json     config_j;
	config_stream >> config_j;
	config_stream.close();
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int num_levels;
	try {
		num_levels = config_j.at("max_levels");
	} catch (nlohmann::detail::out_of_range oor) {
		num_levels = 0;
	}
	double patches_per_proc;
	try {
		patches_per_proc = config_j.at("patches_per_proc");
	} catch (nlohmann::detail::out_of_range oor) {
		patches_per_proc = 0;
	}
	if (num_levels <= 0 || num_levels > t.num_levels) { num_levels = t.num_levels; }
	// generate and balance domain collections
	vector<shared_ptr<SchurHelper<3>>> helpers(num_levels);
	helpers[0] = sh;
	for (int i = 1; i < num_levels; i++) {
		if ((dcs[i]->getGlobalNumDomains() + 0.0) / size < patches_per_proc) {
			num_levels = i;
			break;
		}
		if (dcs[0]->neumann) { dcs[i]->setNeumann(); }
		helpers[i].reset(
		new SchurHelper<3>(*dcs[i], sh->getSolver(), sh->getOp(), sh->getInterpolator()));
	}

	// generate operators
	string op_type;
	try {
		op_type = config_j.at("op_type");
	} catch (nlohmann::detail::out_of_range oor) {
		op_type = "crs_matrix";
	}
	vector<shared_ptr<Operator>> ops(num_levels);
	for (int i = 0; i < num_levels; i++) {
		if (op_type == "crs_matrix") {
			MatrixHelper mh(*dcs[i]);
			ops[i].reset(new MatOp(mh.formCRSMatrix()));
		} else if (op_type == "matrix_free") {
			ops[i].reset(new WrapOp<3>(helpers[i]));
		}
	}

	// generate smoothers
	vector<shared_ptr<Smoother>> smoothers(num_levels);
	smoothers[0].reset(new FFTBlockJacobiSmoother<3>(sh));
	for (int i = 1; i < num_levels; i++) {
		smoothers[i].reset(new FFTBlockJacobiSmoother<3>(helpers[i]));
	}

	// generate inter-level comms, restrictors, interpolators
	vector<shared_ptr<InterLevelComm>> comms(num_levels - 1);
	vector<shared_ptr<Restrictor>>     restrictors(num_levels - 1);
	vector<shared_ptr<Interpolator>>   interpolators(num_levels - 1);
	for (int i = 0; i < num_levels - 1; i++) {
		comms[i].reset(new InterLevelComm(dcs[i + 1], dcs[i]));
		restrictors[i].reset(new AvgRstr(dcs[i + 1], dcs[i], comms[i]));
	}

	// create  level objects
	vector<shared_ptr<Level>> levels(num_levels);
	for (int i = 0; i < num_levels; i++) {
		std::shared_ptr<DCVG> vg(new DCVG(dcs[i]));
		levels[i].reset(new Level(vg));
		levels[i]->setOperator(ops[i]);
		levels[i]->setSmoother(smoothers[i]);
	}

	// set restrictors and interpolators
	string interpolator;
	try {
		interpolator = config_j.at("interpolator");
	} catch (nlohmann::detail::out_of_range oor) {
		interpolator = "trilinear";
	}
	for (int i = 0; i < num_levels - 1; i++) {
		levels[i]->setRestrictor(restrictors[i]);
	}
	if (interpolator == "constant") {
		for (int i = 0; i < num_levels - 1; i++) {
			levels[i + 1]->setInterpolator(
			shared_ptr<Interpolator>(new DrctIntp(dcs[i + 1], dcs[i], comms[i])));
		}
	} else if (interpolator == "trilinear") {
		for (int i = 0; i < num_levels - 1; i++) {
			levels[i + 1]->setInterpolator(
			shared_ptr<Interpolator>(new TriLinIntp(dcs[i + 1], dcs[i], comms[i])));
		}
	} else {
		// TODO throw error
	}

	// link levels to each other
	levels[0]->setCoarser(levels[1]);
	for (int i = 1; i < num_levels - 1; i++) {
		levels[i]->setCoarser(levels[i + 1]);
		levels[i]->setFiner(levels[i - 1]);
	}
	levels[num_levels - 1]->setFiner(levels[num_levels - 2]);

	string cycle_type;
	try {
		cycle_type = config_j.at("cycle_type");
	} catch (nlohmann::detail::out_of_range oor) {
		cycle_type = "V";
	}

	if (cycle_type == "V") {
		cycle.reset(new VCycle(levels[0], config_j));
	} else if (cycle_type == "W") {
		cycle.reset(new WCycle(levels[0], config_j));
	} else {
		// TODO throw error
	}
}
void Helper::apply(Vec f, Vec u)
{
	cycle->apply(f, u);
}
