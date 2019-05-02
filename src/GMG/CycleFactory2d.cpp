#include "CycleFactory2d.h"
#include <GMG/AvgRstr.h>
#include <GMG/DrctIntp.h>
#include <GMG/FFTBlockJacobiSmoother.h>
#include <GMG/InterLevelComm.h>
#include <GMG/VCycle.h>
#include <GMG/WCycle.h>
#include <MatrixHelper2d.h>
#include <Operators/PetscMatOp.h>
#include <Operators/SchurDomainOp.h>
#include <json.hpp>
#include <limits>
using namespace GMG;
using namespace std;
using nlohmann::json;
static std::shared_ptr<Operator<2>> getNewOperator(std::string                          op_type,
                                                   std::shared_ptr<DomainCollection<2>> dc,
                                                   std::shared_ptr<SchurHelper<2>>      sh)
{
	std::shared_ptr<Operator<2>> op;
	if (op_type == "crs_matrix") {
		MatrixHelper2d mh(*dc);
		op.reset(new PetscMatOp<2>(mh.formCRSMatrix()));
	} else if (op_type == "matrix_free") {
		op.reset(new SchurDomainOp<2>(sh));
	}
	return op;
}
static std::shared_ptr<Smoother<2>> getNewSmoother(std::string smoother_type,
                                                   std::shared_ptr<DomainCollection<2>> dc,
                                                   std::shared_ptr<SchurHelper<2>>      sh)
{
	return std::shared_ptr<Smoother<2>>(new FFTBlockJacobiSmoother<2>(sh));
}
static std::shared_ptr<Interpolator<2>>
getNewInterpolator(std::string interpolator_type, std::shared_ptr<DomainCollection<2>> dc,
                   std::shared_ptr<DomainCollection<2>> finer_dc,
                   std::shared_ptr<InterLevelComm<2>>   ilc)
{
	return std::shared_ptr<Interpolator<2>>(new DrctIntp<2>(dc, finer_dc, ilc));
}
static std::shared_ptr<Restrictor<2>>
getNewRestrictor(std::string restrictor_type, std::shared_ptr<DomainCollection<2>> dc,
                 std::shared_ptr<DomainCollection<2>> finer_dc,
                 std::shared_ptr<InterLevelComm<2>>   ilc)
{
	return std::shared_ptr<Restrictor<2>>(new AvgRstr<2>(dc, finer_dc, ilc));
}
std::shared_ptr<Cycle<2>>
CycleFactory2d::getCycle(std::string config_file, std::shared_ptr<DomainCollectionGenerator<2>> dcg,
                         std::shared_ptr<PatchSolver<2>>   solver,
                         std::shared_ptr<PatchOperator<2>> op,
                         std::shared_ptr<IfaceInterp<2>>   interp)
{
	ifstream config_stream(config_file);
	json     config_j;
	config_stream >> config_j;
	config_stream.close();
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int max_levels;
	try {
		max_levels = config_j.at("max_levels");
	} catch (nlohmann::detail::out_of_range oor) {
		max_levels = 0;
	}
	double patches_per_proc;
	try {
		patches_per_proc = config_j.at("patches_per_proc");
	} catch (nlohmann::detail::out_of_range oor) {
		patches_per_proc = 0;
	}

	string op_type;
	try {
		op_type = config_j.at("op_type");
	} catch (nlohmann::detail::out_of_range oor) {
		op_type = "crs_matrix";
	}

	string smoother_type = "block jacobi";
	// set restrictors and interpolators
	string interpolator_type;
	try {
		interpolator_type = config_j.at("interpolator");
	} catch (nlohmann::detail::out_of_range oor) {
		interpolator_type = "trilinear";
	}
	string restrictor_type;

	// finest level
	shared_ptr<Level<2>>            finest_level;
	shared_ptr<DomainCollection<2>> finer_dc;
	{
		shared_ptr<DomainCollection<2>> dc = dcg->getFinestDC();
		shared_ptr<SchurHelper<2>>      sh(new SchurHelper<2>(dc, solver, op, interp));
		shared_ptr<VectorGenerator<2>>  vg(new DomainCollectionVG<2>(dc));
		finest_level.reset(new Level<2>(vg));
		finest_level->setOperator(getNewOperator(op_type, dc, sh));
		finest_level->setSmoother(getNewSmoother(smoother_type, dc, sh));

		finer_dc = dc;
	}
	shared_ptr<Level<2>> finer_level = finest_level;
	// other levels
	int curr_level = 1;
	while (dcg->hasCoarserDC() && (max_levels <= 0 || curr_level < max_levels)) {
		// create new level
		shared_ptr<DomainCollection<2>> dc = dcg->getCoarserDC();
		if ((dc->getGlobalNumDomains() + 0.0) / size < patches_per_proc) { break; }
		shared_ptr<SchurHelper<2>>     sh(new SchurHelper<2>(dc, solver, op, interp));
		shared_ptr<VectorGenerator<2>> vg(new DomainCollectionVG<2>(dc));
		shared_ptr<Level<2>>           coarser_level(new Level<2>(vg));
		coarser_level->setOperator(getNewOperator(op_type, dc, sh));
		coarser_level->setSmoother(getNewSmoother(smoother_type, dc, sh));

		// link levels
		coarser_level->setFiner(finer_level);
		finer_level->setCoarser(coarser_level);

		// set restrictor and interpolator operators
		shared_ptr<InterLevelComm<2>> ilc(new InterLevelComm<2>(dc, finer_dc));
		finer_level->setRestrictor(getNewRestrictor(restrictor_type, dc, finer_dc, ilc));
		coarser_level->setInterpolator(getNewInterpolator(interpolator_type, dc, finer_dc, ilc));

		curr_level++;
		finer_level = coarser_level;
		finer_dc    = dc;
	}
	string cycle_type;
	try {
		cycle_type = config_j.at("cycle_type");
	} catch (nlohmann::detail::out_of_range oor) {
		cycle_type = "V";
	}
	shared_ptr<Cycle<2>> cycle;
	if (cycle_type == "V") {
		cycle.reset(new VCycle<2>(finest_level, config_j));
	} else if (cycle_type == "W") {
		cycle.reset(new WCycle<2>(finest_level, config_j));
	} else {
		throw 3;
	}
	return cycle;
}
