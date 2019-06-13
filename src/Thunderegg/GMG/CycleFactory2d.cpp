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

#include "CycleFactory2d.h"
#include <Thunderegg/GMG/AvgRstr.h>
#include <Thunderegg/GMG/DrctIntp.h>
#include <Thunderegg/GMG/FFTBlockJacobiSmoother.h>
#include <Thunderegg/GMG/InterLevelComm.h>
#include <Thunderegg/GMG/VCycle.h>
#include <Thunderegg/GMG/WCycle.h>
#include <Thunderegg/MatrixHelper2d.h>
#include <Thunderegg/Operators/PetscMatOp.h>
#include <Thunderegg/Operators/SchurDomainOp.h>
#include <limits>
#include <locale>
using namespace GMG;
using namespace std;
static std::shared_ptr<Operator<2>> getNewOperator(std::string                     op_type,
                                                   std::shared_ptr<Domain<2>>      domain,
                                                   std::shared_ptr<SchurHelper<2>> sh)
{
	std::shared_ptr<Operator<2>> op;
	if (op_type == "crs_matrix") {
		MatrixHelper2d mh(domain);
		op.reset(new PetscMatOp<2>(mh.formCRSMatrix()));
	} else if (op_type == "matrix_free") {
		op.reset(new SchurDomainOp<2>(sh));
	}
	return op;
}
static std::shared_ptr<Smoother<2>> getNewSmoother(std::string                     smoother_type,
                                                   std::shared_ptr<Domain<2>>      domain,
                                                   std::shared_ptr<SchurHelper<2>> sh)
{
	return std::shared_ptr<Smoother<2>>(new FFTBlockJacobiSmoother<2>(sh));
}
static std::shared_ptr<Interpolator<2>> getNewInterpolator(std::string interpolator_type,
                                                           std::shared_ptr<Domain<2>> domain,
                                                           std::shared_ptr<Domain<2>> finer_domain,
                                                           std::shared_ptr<InterLevelComm<2>> ilc)
{
	return std::shared_ptr<Interpolator<2>>(new DrctIntp<2>(domain, finer_domain, ilc));
}
static std::shared_ptr<Restrictor<2>> getNewRestrictor(std::string                restrictor_type,
                                                       std::shared_ptr<Domain<2>> domain,
                                                       std::shared_ptr<Domain<2>> finer_domain,
                                                       std::shared_ptr<InterLevelComm<2>> ilc)
{
	return std::shared_ptr<Restrictor<2>>(new AvgRstr<2>(domain, finer_domain, ilc));
}
std::shared_ptr<Cycle<2>>
CycleFactory2d::getCycle(const CycleOpts &opts, std::shared_ptr<DomainCollectionGenerator<2>> dcg,
                         std::shared_ptr<PatchSolver<2>>   solver,
                         std::shared_ptr<PatchOperator<2>> op,
                         std::shared_ptr<IfaceInterp<2>>   interp)
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	string op_type = "matrix_free";

	string smoother_type = "block jacobi";
	// set restrictors and interpolators
	string interpolator_type = "constant";
	string restrictor_type;

	// finest level
	shared_ptr<Level<2>>  finest_level;
	shared_ptr<Domain<2>> finer_domain;
	{
		shared_ptr<Domain<2>>          domain = dcg->getFinestDC();
		shared_ptr<SchurHelper<2>>     sh(new SchurHelper<2>(domain, solver, op, interp));
		shared_ptr<VectorGenerator<2>> vg(new DomainVG<2>(domain));
		finest_level.reset(new Level<2>(vg));
		finest_level->setOperator(getNewOperator(op_type, domain, sh));
		finest_level->setSmoother(getNewSmoother(smoother_type, domain, sh));

		finer_domain = domain;
	}
	shared_ptr<Level<2>> finer_level = finest_level;
	// other levels
	int curr_level = 1;
	while (dcg->hasCoarserDC() && (opts.max_levels <= 0 || curr_level < opts.max_levels)) {
		// create new level
		shared_ptr<Domain<2>> domain = dcg->getCoarserDC();
		if ((domain->getNumGlobalPatches() + 0.0) / size < opts.patches_per_proc) { break; }
		shared_ptr<SchurHelper<2>>     sh(new SchurHelper<2>(domain, solver, op, interp));
		shared_ptr<VectorGenerator<2>> vg(new DomainVG<2>(domain));
		shared_ptr<Level<2>>           coarser_level(new Level<2>(vg));
		coarser_level->setOperator(getNewOperator(op_type, domain, sh));
		coarser_level->setSmoother(getNewSmoother(smoother_type, domain, sh));

		// link levels
		coarser_level->setFiner(finer_level);
		finer_level->setCoarser(coarser_level);

		// set restrictor and interpolator operators
		shared_ptr<InterLevelComm<2>> ilc(new InterLevelComm<2>(domain, finer_domain));
		finer_level->setRestrictor(getNewRestrictor(restrictor_type, domain, finer_domain, ilc));
		coarser_level->setInterpolator(
		getNewInterpolator(interpolator_type, domain, finer_domain, ilc));

		curr_level++;
		finer_level  = coarser_level;
		finer_domain = domain;
	}
	shared_ptr<Cycle<2>> cycle;
	if (opts.cycle_type == "V") {
		cycle.reset(new VCycle<2>(finest_level, opts));
	} else if (opts.cycle_type == "W") {
		cycle.reset(new WCycle<2>(finest_level, opts));
	} else {
		throw 3;
	}
	return cycle;
}
