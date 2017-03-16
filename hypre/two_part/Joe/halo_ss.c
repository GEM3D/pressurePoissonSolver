





#include <math.h>
#include <stdio.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_parcsr_ls.h"
#include "vis.c"
#include "add_to_vis.c"
#define PI2 M_PI *M_PI
#include "HYPRE_parcsr_mv.h"
#include "mpi.h"

int main(int argc, char *argv[]) {

	int i,j;
	int myid, num_procs;
	int n;
	//int vis;
	int ilower[2], iupper[2];
	int nparts = 2;
	int part0 = 0;
	int part1 = 1;
	int nvars = 1;
	int var  = 0;
	int num_stencil_entries = 5;
	int entry;
	int nvalues;
	double *values;
	int stencil_indices[5];
	int object_type = HYPRE_PARCSR;
	int num_iterations;
	double final_res_norm, t1, t2;
	double h;
	double h2inv;	
	double X, Y;

printf("variables setup\n");

	HYPRE_SStructGraph graph;
  	HYPRE_SStructGrid grid;
  	HYPRE_SStructStencil stencil;
  	HYPRE_SStructMatrix A;
  	HYPRE_ParCSRMatrix parcsr_A;
 	HYPRE_SStructVector b;
  	HYPRE_ParVector par_b;
  	HYPRE_SStructVector x;
  	HYPRE_ParVector par_x;

	HYPRE_Solver solver;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


	//Set defaults
	n = 100;
	//vis = 1;
	h = 1.0/(n);
	h2inv = 1.0/(h*h);


printf("defaults set\n");

	ilower[0] = 1;
	ilower[1] = 1;

	iupper[0] = ilower[0] + n-1;
	iupper[1] = ilower[1] + n-1;


	HYPRE_SStructGridCreate(MPI_COMM_WORLD, 2, nparts, &grid);
	
	HYPRE_SStructGridSetExtents(grid, part0, ilower, iupper);

printf("first\n");


	HYPRE_SStructGridSetExtents(grid, part1, ilower, iupper);

printf("grid extents set\n");

	HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
	
	HYPRE_SStructGridSetVariables(grid, part0, nvars, vartypes);
	HYPRE_SStructGridSetVariables(grid, part1, nvars, vartypes);

//Set Neighbor Part
	int tupper[2] = {n,n};
	int tlower[2] = {n,1};
	int imap[2] = {0,1};
	int idir[2] = {1,1};
	HYPRE_SStructGridSetNeighborPart(grid, part0, tlower, tupper, part1, tlower, tupper, imap, idir);
	

	tupper[0] = 1;
	tupper[1] = n;
	tlower[0] = 1;
	tlower[1] = 1;
	HYPRE_SStructGridSetNeighborPart(grid, part1, tlower, tupper, part0, tlower, tupper, imap, idir);

	
	HYPRE_SStructGridAssemble(grid);

printf("grid set\n");

	HYPRE_SStructStencilCreate(2, 5, &stencil);
	int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0, -1}, {0, 1}};

	
	for (entry = 0; entry < 5; entry++) {
		HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], var);
	}
	

printf("Why\n");
	/*
	HYPRE_SStructStencilSetEntry(stencil, 0, offsets[0], var);
	HYPRE_SStructStencilSetEntry(stencil, 1, offsets[1], var);
	HYPRE_SStructStencilSetEntry(stencil, 2, offsets[2], var);
	HYPRE_SStructStencilSetEntry(stencil, 3, offsets[3], var);
	HYPRE_SStructStencilSetEntry(stencil, 4, offsets[4], var);
	*/
printf("What\n");
	HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);

printf("Did\n");

	HYPRE_SStructGraphSetObjectType(graph, object_type);
	
printf("I\n");
	HYPRE_SStructGraphSetStencil(graph, part1, var, stencil);

printf("To\n");

	//HYPRE_SStructGraphSetStencil(graph, part1, var, stencil);
	HYPRE_SStructGraphSetStencil(graph, part0, var, stencil);

printf("Deserve\n");

printf("halo begin\n");
	
	//Halo Stuff

        HYPRE_Int *halo;
        halo = calloc(n, sizeof(int));
        //for(j = 0; j < n; j++){
	//	halo[0] = j*99+99;
	//	HYPRE_SStructGraphAddEntries(graph, 1, halo, var, 0, halo, var);
        //}
	
        for(j = 0; j < n; j+=1) {
                halo[j] = j*99+99;
        }
printf("forloop\n");
        HYPRE_SStructGraphAddEntries(graph, 1, halo, var, 0, halo, var);
        HYPRE_SStructGraphAddEntries(graph, 0, halo, var, 1, halo, var);
	



	HYPRE_SStructGraphAssemble(graph);

printf("graph set\n");

	nvalues = num_stencil_entries * n * n;
	
	HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
	HYPRE_SStructMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_SStructMatrixInitialize(A);
		
	values = calloc(nvalues, sizeof(double));
	for (i = 0; i < num_stencil_entries; i++);{
		stencil_indices[i] = i;
		printf("%d\n",stencil_indices[i]);
	}
	

	stencil_indices[0] = 0;
	stencil_indices[1] = 1;
	stencil_indices[2] = 2;
	stencil_indices[3] = 3;
	stencil_indices[4] = 4;

printf("stencil indices set\n");

	values = calloc(nvalues, sizeof(double));
	for (i = 0; i < nvalues; i += num_stencil_entries) {
		values[i] = -4.0 * h2inv;
		values[i+1] = h2inv;
		values[i+2] = h2inv;
		values[i+3] = h2inv;
		values[i+4] = h2inv;
	}

printf("matrix values pre-set\n");
printf("%d\n", nvalues);
printf("[%d,%d] [%d, %d]\n", ilower[0], ilower[1], iupper[0], iupper[1]);
printf("part: %d, var: %d, numstencilentries: %d\n", part0, var, num_stencil_entries);
	HYPRE_SStructMatrixSetBoxValues(A, part0, ilower, iupper, var, num_stencil_entries, stencil_indices, values);

printf("matrix values set part0\n");

	HYPRE_SStructMatrixSetBoxValues(A, part1, ilower, iupper, var, num_stencil_entries, stencil_indices, values);

printf("matrix set\n");

	//Set Neumann boundary conditions on both parts. Assumer part0 lies to the left of part1;
	
	int bc_ilower[2];
	int bc_iupper[2];
	int nentries = 1;
	int nvaluesx = nentries * n;
	int nvaluesy = nentries * n;
	double *valuesx, *valuesy, *center_valuesx, *center_valuesy;
	int stencil_index[1];
	int center_index[1];
	center_index[0] = 0;

	valuesx = calloc(nvaluesx, sizeof(double));
	valuesy = calloc(nvaluesy, sizeof(double));
	center_valuesx = calloc(nvaluesx, sizeof(double));
	center_valuesy = calloc(nvaluesy, sizeof(double));
	for (i = 0; i< nvaluesx; i++){
		valuesx[i] = 0.0;
		valuesy[i] = 0.0;
		center_valuesx[i] = h2inv;
		center_valuesy[i] = h2inv;
	}

	//Bottom row of gridpoints
	bc_ilower[0] = 1;
	bc_ilower[1] = 1;
	
	bc_iupper[0] = n;
	bc_iupper[1] = 1;

	stencil_index[0] = 3;

	center_valuesy[3] = -h2inv;
	
	HYPRE_SStructMatrixSetBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesy);
	HYPRE_SStructMatrixAddToBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, center_index, valuesy);
	center_valuesy[3] = h2inv;

	HYPRE_SStructMatrixSetBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesy); 
	HYPRE_SStructMatrixAddToBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, center_index, valuesy);

	//Top row of gridpoints
	bc_ilower[0] = 1;
	bc_ilower[1] = n;

	bc_iupper[0] = n;
	bc_iupper[1] = n;

	stencil_index[0] = 4;
	
	HYPRE_SStructMatrixSetBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
	HYPRE_SStructMatrixAddToBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, center_index, valuesx);
	HYPRE_SStructMatrixSetBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
	HYPRE_SStructMatrixAddToBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, center_index, valuesx);

	//Left Column of gridpoints (Note only gets set for part0 as this portion of part1 is not on the boundary)
	bc_ilower[0] = 1;
	bc_ilower[1] = 1;
	
	bc_iupper[0] = 1;
	bc_iupper[1] = n;

	stencil_index[0] = 1;
	
	HYPRE_SStructMatrixSetBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
	HYPRE_SStructMatrixAddToBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, center_index, valuesx);
	
	//Right column of gridpoints (Note only gets set for part1 as this portion of part0 is not on the boundary)
	bc_ilower[0] = n;
	bc_ilower[1] = 1;
	
	bc_iupper[0] = n;
	bc_iupper[1] = n;

	stencil_index[0] = 2;
	
	HYPRE_SStructMatrixSetBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);	
	HYPRE_SStructMatrixAddToBoxValues(A, part1,bc_ilower, bc_iupper, var, nentries, center_index, valuesx);
	free(valuesx);
	free(center_valuesx);
	free(valuesy);
	free(center_valuesy);

	HYPRE_SStructMatrixAssemble(A);
	HYPRE_SStructMatrixGetObject(A, (void **)&parcsr_A);

printf("boundary conditions set\n");

	//Create the right hand side and solution vectors
	HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &b);
	HYPRE_SStructVectorSetObjectType(b, HYPRE_PARCSR);
	HYPRE_SStructVectorInitialize(b);

	HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &x);
	HYPRE_SStructVectorSetObjectType(x, HYPRE_PARCSR);
	HYPRE_SStructVectorInitialize(x);

	//Set the right hand side, exact solution, and initial guess values
	double *rhs_values, *x_values, *exactsolution;
	nvalues = n * n;

	rhs_values = calloc(nvalues, sizeof(double));
	x_values = calloc(nvalues, sizeof(double));
	exactsolution = calloc(nvalues, sizeof(double));

	for (j = 0; j < n; j++) {
		Y = (j + .5) * h;
		for(i = 0; i < n; i++) {
			X = (i + .5) * h;
			rhs_values[i + j * n] = -8. * PI2 * cos(2. * M_PI * X)*cos(2 * M_PI * Y);
			exactsolution[i + j * n] = 1.0 * cos(2. * M_PI * X) * cos(2. * M_PI * Y);
			x_values[i+ j * n] = 0.0;
		}
	}

	HYPRE_SStructVectorSetBoxValues(b, part0, ilower, iupper, var, rhs_values);
	HYPRE_SStructVectorSetBoxValues(x, part0, ilower, iupper, var, x_values);
	HYPRE_SStructVectorSetBoxValues(x, part1, ilower, iupper, var, x_values);

        for (j = 0; j < n; j++) {
               Y = (j + .5) * h;
               for(i = 0; i < n; i++) {
                       X = (n*part1 * i + .5) * h;
                       rhs_values[i + j * n] = -8. * PI2 * cos(2. * M_PI * X)*cos(2 * M_PI * Y);
                       exactsolution[i + j * n] = 1.0 * cos(2. * M_PI * X) * cos(2. * M_PI * Y);
                       x_values[i+ j * n] = 0.0;
                }
        }

	HYPRE_SStructVectorSetBoxValues(b, part1, ilower, iupper, var, rhs_values);

printf("rhs and x vectors set\n");

	free(x_values);
	free(rhs_values);

	HYPRE_SStructVectorAssemble(b);
	HYPRE_SStructVectorGetObject(b, (void **)&par_b);
	HYPRE_SStructVectorAssemble(x);
	HYPRE_SStructVectorGetObject(x, (void **)&par_x);

	HYPRE_SStructMatrixPrint("ss_halo_data/ss.initial.A", A, 0);
	HYPRE_SStructVectorPrint("ss_halo_data/ss.initial.A", b, 0);


	HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
	HYPRE_ParCSRPCGSetTol(solver, 1.0e-10);
	HYPRE_ParCSRPCGSetMaxIter(solver, 500);
	HYPRE_ParCSRPCGSetTwoNorm(solver, 1);
	HYPRE_ParCSRPCGSetLogging(solver, 3);
	
	HYPRE_Solver precond;

	HYPRE_BoomerAMGCreate(&precond);
	HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
	HYPRE_BoomerAMGSetCoarsenType(precond, 6);
        HYPRE_BoomerAMGSetTol(precond, 1.0e-12);
        HYPRE_BoomerAMGSetPrintLevel(precond, 1);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0);
        HYPRE_BoomerAMGSetMaxIter(precond, 5);
        HYPRE_BoomerAMGSetNumSweeps(precond, 2);
        HYPRE_BoomerAMGSetRelaxType(precond, 0);
        HYPRE_BoomerAMGSetRelaxWt(precond, .3);


	HYPRE_ParCSRPCGSetPrecond(solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
	HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
	t1 = MPI_Wtime();
	HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);
	t2 = MPI_Wtime();

	double time = t2-t1;

	HYPRE_ParCSRPCGGetNumIterations(solver, &num_iterations);
	HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
	
	printf("\n");
	printf("Iterations = %d\n", num_iterations);
	printf("Final Relative Residual Norm = %e\n", final_res_norm);
	printf("Elapsed Wall Time = %f seconds\n", time);
	printf("\n");
	

	char solutionfile[255];
	sprintf(solutionfile, "%s.%06d", "ss_halo_data.final.x", myid);
	HYPRE_ParVectorPrint(par_x, solutionfile);

	FILE *solution;
	sprintf(solutionfile, "%s.%06d.%d", "ss_halo_data.final.x", myid, myid);
	free(values);
	double totalvalues = nparts * n * n;
	double diff = 0.0;
	double diff2;
	double sum = 0.0;
	values = calloc(totalvalues+1, sizeof(double));


	if ((solution = fopen(solutionfile, "rt")) == NULL) {
		printf("Error: can't open output file %s\n", solutionfile);
		MPI_Finalize();
		exit(1);
	}
	double mean = 0.0;

	for (i = 0; i < totalvalues + 1; i++) {
		fscanf(solution, "%lf", &values[i]);
		if (i != 0){
			mean += values[i];
		}		
	}

	mean /= nvalues;
	printf("The mean is: %f\n", mean);
	for (i = 1; i < totalvalues + 1; i++) {
		diff = values[i] -mean - exactsolution[i-1];
		diff2 = diff * diff;
		if (diff2 > 10){
			printf("%+6f       %+6f     %+6f\n", values[i], exactsolution[i-1],diff);
		}
		sum += diff2;
	}	

	double L2 = sum/totalvalues;

	printf("The L2 norm is %e\n", L2);
	
	fflush(solution);
	fclose(solution);


	HYPRE_BoomerAMGDestroy(precond);
	HYPRE_ParCSRPCGDestroy(solver);


	HYPRE_SStructMatrixDestroy(A);
	HYPRE_SStructVectorDestroy(b);
	HYPRE_SStructVectorDestroy(x);

	MPI_Finalize();

	return(0);
}




