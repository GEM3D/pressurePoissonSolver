/***************************************************************************
 *  Thunderegg, a library for solving Poisson's equation on adaptively
 *  refined block-structured Cartesian grids
 *
 *  Copyright (C) 2019  Thunderegg Developers. See AUTHORS.md file at the
 *  top-level directory.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ***************************************************************************/

#include "CycleFactory3d.h"
#include <Thunderegg/GMG/AvgRstr.h>
#include <Thunderegg/GMG/DrctIntp.h>
#include <Thunderegg/GMG/FFTBlockJacobiSmoother.h>
#include <Thunderegg/GMG/InterLevelComm.h>
#include <Thunderegg/GMG/VCycle.h>
#include <Thunderegg/GMG/WCycle.h>
#include <Thunderegg/MatrixHelper.h>
#include <Thunderegg/Operators/PetscMatOp.h>
#include <Thunderegg/Operators/SchurDomainOp.h>
#include <fstream>
#include <limits>
using namespace GMG;
using namespace std;
static std::shared_ptr<Operator<3>> getNewOperator(std::string                     op_type,
                                                   std::shared_ptr<Domain<3>>      domain,
                                                   std::shared_ptr<SchurHelper<3>> sh)
{
	std::shared_ptr<Operator<3>> op;
	if (op_type == "crs_matrix") {
		MatrixHelper mh(domain);
		op.reset(new PetscMatOp<3>(mh.formCRSMatrix()));
	} else if (op_type == "matrix_free") {
		op.reset(new SchurDomainOp<3>(sh));
	}
	return op;
}
static std::shared_ptr<Smoother<3>> getNewSmoother(std::string                     smoother_type,
                                                   std::shared_ptr<Domain<3>>      domain,
                                                   std::shared_ptr<SchurHelper<3>> sh)
{
	return std::shared_ptr<Smoother<3>>(new FFTBlockJacobiSmoother<3>(sh));
}
static std::shared_ptr<Interpolator<3>> getNewInterpolator(std::string interpolator_type,
                                                           std::shared_ptr<Domain<3>> domain,
                                                           std::shared_ptr<Domain<3>> finer_domain,
                                                           std::shared_ptr<InterLevelComm<3>> ilc)
{
	return std::shared_ptr<Interpolator<3>>(new DrctIntp<3>(domain, finer_domain, ilc));
}
static std::shared_ptr<Restrictor<3>> getNewRestrictor(std::string                restrictor_type,
                                                       std::shared_ptr<Domain<3>> domain,
                                                       std::shared_ptr<Domain<3>> finer_domain,
                                                       std::shared_ptr<InterLevelComm<3>> ilc)
{
	return std::shared_ptr<Restrictor<3>>(new AvgRstr<3>(domain, finer_domain, ilc));
}
std::shared_ptr<Cycle<3>> CycleFactory3d::getCycle(const CycleOpts &                   opts,
                                                   std::shared_ptr<DomainGenerator<3>> dcg,
                                                   std::shared_ptr<PatchSolver<3>>     solver,
                                                   std::shared_ptr<PatchOperator<3>>   op,
                                                   std::shared_ptr<IfaceInterp<3>>     interp)
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	string op_type = "matrix_free";

	string smoother_type = "block jacobi";
	// set restrictors and interpolators
	string interpolator_type = "constant";
	string restrictor_type;

	// finest level
	shared_ptr<Level<3>>  finest_level;
	shared_ptr<Domain<3>> finer_domain;
	{
		shared_ptr<Domain<3>>          domain = dcg->getFinestDomain();
		shared_ptr<SchurHelper<3>>     sh(new SchurHelper<3>(domain, solver, op, interp));
		shared_ptr<VectorGenerator<3>> vg(new DomainVG<3>(domain));
		finest_level.reset(new Level<3>(vg));
		finest_level->setOperator(getNewOperator(op_type, domain, sh));
		finest_level->setSmoother(getNewSmoother(smoother_type, domain, sh));

		finer_domain = domain;
	}
	shared_ptr<Level<3>> finer_level = finest_level;
	// other levels
	int curr_level = 1;
	while (dcg->hasCoarserDomain() && (opts.max_levels <= 0 || curr_level < opts.max_levels)) {
		// create new level
		shared_ptr<Domain<3>> domain = dcg->getCoarserDomain();
		if ((domain->getNumGlobalPatches() + 0.0) / size < opts.patches_per_proc) { break; }
		shared_ptr<SchurHelper<3>>     sh(new SchurHelper<3>(domain, solver, op, interp));
		shared_ptr<VectorGenerator<3>> vg(new DomainVG<3>(domain));
		shared_ptr<Level<3>>           coarser_level(new Level<3>(vg));
		coarser_level->setOperator(getNewOperator(op_type, domain, sh));
		coarser_level->setSmoother(getNewSmoother(smoother_type, domain, sh));

		// link levels
		coarser_level->setFiner(finer_level);
		finer_level->setCoarser(coarser_level);

		// set restrictor and interpolator operators
		shared_ptr<InterLevelComm<3>> ilc(new InterLevelComm<3>(domain, finer_domain));
		finer_level->setRestrictor(getNewRestrictor(restrictor_type, domain, finer_domain, ilc));
		coarser_level->setInterpolator(
		getNewInterpolator(interpolator_type, domain, finer_domain, ilc));

		curr_level++;
		finer_level  = coarser_level;
		finer_domain = domain;
	}
	shared_ptr<Cycle<3>> cycle;
	if (opts.cycle_type == "V") {
		cycle.reset(new VCycle<3>(finest_level, opts));
	} else if (opts.cycle_type == "W") {
		cycle.reset(new WCycle<3>(finest_level, opts));
	} else {
		throw 3;
	}
	return cycle;
}
