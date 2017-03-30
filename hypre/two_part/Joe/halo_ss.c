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

int main(int argc, char *argv[]) 
{

	int i ,j;
	int myid, num_procs;
	int n;
	int ilower[2], iupper[2];
	int nparts = 2;
	int part0  = 0;
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
    HYPRE_Solver   precond;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


	//Set defaults
	n = 64;
	h = 1.0/(n);
	h2inv = 1.0/(h*h);


	printf("defaults set\n");

	ilower[0] = 1;
	ilower[1] = 1;

	iupper[0] = n;
	iupper[1] = n;

	/* Added by Donna */
	int ll[2], lr[2], ul[2], ur[2];
	ll[0] = 1;
	ll[1] = 1;
	lr[0] = n;
	lr[1] = 1;
	ul[0] = 1;
	ul[1] = n;
	ur[0] = n;
	ur[1] = n;




	HYPRE_SStructGridCreate(MPI_COMM_WORLD, 2, nparts, &grid);
	
	HYPRE_SStructGridSetExtents(grid, part0, ilower, iupper);

	printf("first\n");


	HYPRE_SStructGridSetExtents(grid, part1, ilower, iupper);

	printf("grid extents set\n");

	HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
	
	HYPRE_SStructGridSetVariables(grid, part0, nvars, vartypes);
	HYPRE_SStructGridSetVariables(grid, part1, nvars, vartypes);



	//temporary set neighbor part for testing
	//part 0 to 1
	{
		int lr_g[2], ur_g[2], ll_g[2], ul_g[2];
		int part, nbor_part;
		int index_map[2] = {0,1};
		int index_dir[2] = {1,1};

		lr_g[0] = lr[0] + 1;
		lr_g[1] = lr[1];
		ur_g[0] = ur[0] + 1;
		ur_g[1] = ur[1];

 #if 0	
	int b_ilower[2] = {n,1}; 
	int b_iupper[2] = {n,n};
	int nbor_ilower[2] = {1,1}; 
	int nbor_iupper[2] = {1,n};
#endif	


		part = 0;
		nbor_part = 1;
		HYPRE_SStructGridSetNeighborPart(grid, part, lr_g, ur_g, nbor_part, ll, ul,
		                                 index_map, index_dir);


		ll_g[0] = ll[0] - 1;
		ll_g[1] = ll[1];
		ul_g[0] = ul[0] - 1;
		ul_g[1] = ul[1];


		part = 1;
		nbor_part = 0;
		HYPRE_SStructGridSetNeighborPart(grid, part, ll_g, ul_g, part, lr, ur, index_map, index_dir);

		HYPRE_SStructGridAssemble(grid);

	}

	printf("grid set\n");

	/* ------------------------------------------------- 
		Create stencil
		------------------------------------------------- */
	{
		HYPRE_SStructStencilCreate(2, 5, &stencil);

		int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0, -1}, {0, 1}};


		for (entry = 0; entry < 5; entry++) 
		{	
			HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], var);
		}	
	}

    /* -------------------------------------------------------------
       3. Set up the graph.
       ------------------------------------------------------------- */
	{
		HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);
		HYPRE_SStructGraphSetObjectType(graph, object_type);

		for(int part = 0; part < nparts; part++)
		{
			HYPRE_SStructGraphSetStencil(graph, part, var, stencil);
		}

		printf("halo begin\n");


		//Halo Stuff
#if 0
		int *halo1;
		int *halo2;
		halo1 = calloc(2, sizeof(int));
		halo2 = calloc(2, sizeof(int));

		for(j = 1; j <= n; j++) 
		{
			halo1[0] = n- 1 ;
			halo1[1] = j;
			halo2[0] = 2;
			halo2[1] = j;
			//printf("%d, %d\n", halo1[1], halo2[1]);

			HYPRE_SStructGraphAddEntries(graph, 1, halo2, var, 0, halo1, var);
			HYPRE_SStructGraphAddEntries(graph, 0, halo1, var, 1, halo2, var);
		}		
#endif 		
		HYPRE_SStructGraphAssemble(graph);
	}


	printf("graph set\n");

    /* -------------------------------------------------------------
       4. Set up a SStruct Matrix
       ------------------------------------------------------------- */
	{
		nvalues = num_stencil_entries * n * n;

		HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
		HYPRE_SStructMatrixSetObjectType(A, HYPRE_PARCSR);
		HYPRE_SStructMatrixInitialize(A);

		values = calloc(nvalues, sizeof(double));

		for (i = 0; i < num_stencil_entries; i++){
			stencil_indices[i] = i;
		}


		printf("stencil indices set\n");

		values = calloc(nvalues, sizeof(double));
		for (i = 0; i < nvalues; i += num_stencil_entries) 
		{
			values[i] = -4.0 * h2inv;
			values[i+1] = h2inv;
			values[i+2] = h2inv;
			values[i+3] = h2inv;
			values[i+4] = h2inv;	
		}

		printf("matrix values pre-set\n");


		HYPRE_SStructMatrixSetBoxValues(A, part0, ilower, iupper, var, 
		                                num_stencil_entries, stencil_indices, values);

		printf("matrix values set part0\n");

		HYPRE_SStructMatrixSetBoxValues(A, part1, ilower, iupper, var, 
		                                num_stencil_entries, stencil_indices, values);

		printf("matrix set\n");
	}	

	//Set Neumann boundary conditions on both parts. Assume part0 lies to the left of part1;
	
	{
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
		printf("begin nvalues loop\n");

		for (i = 0; i< nvaluesx; i++)
		{
			valuesx[i] = 0.0;
			valuesy[i] = 0.0;
			center_valuesx[i] = h2inv;
			center_valuesy[i] = h2inv;
		}

		printf("bc begin set\n");

		//Bottom row of gridpoints
		bc_ilower[0] = 1;
		bc_ilower[1] = 1;

		bc_iupper[0] = n;
		bc_iupper[1] = 1;

		stencil_index[0] = 3;

		center_valuesy[1] = -h2inv;

		printf("bottom part 0\n");

		HYPRE_SStructMatrixSetBoxValues(A, part0, bc_ilower, bc_iupper, 
		                                var, nentries, stencil_index, valuesy);

		printf("set\n");
		HYPRE_SStructMatrixAddToBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesy);
		center_valuesy[1] = h2inv;

		printf("bottom part 1\n");

		HYPRE_SStructMatrixSetBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesy); 
		//HYPRE_SStructMatrixAddToBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesy);
		HYPRE_SStructMatrixAddToBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesy);	
		center_valuesy[1] = h2inv;
		//printf("bottom row set\n");

		//Top row of gridpoints
		bc_ilower[0] = 1;
		bc_ilower[1] = n;

		bc_iupper[0] = n;
		bc_iupper[1] = n;

		stencil_index[0] = 4;

		HYPRE_SStructMatrixSetBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
		HYPRE_SStructMatrixAddToBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);
		HYPRE_SStructMatrixSetBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
		HYPRE_SStructMatrixAddToBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);

		//printf("top row set\n");

		//Left Column of gridpoints (Note only gets set for part0 as this portion of part1 is not on the boundary)
		bc_ilower[0] = 1;
		bc_ilower[1] = 1;

		bc_iupper[0] = 1;
		bc_iupper[1] = n;

		stencil_index[0] = 1;

		HYPRE_SStructMatrixSetBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
		HYPRE_SStructMatrixAddToBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);

		//REMOVE THIS
		//	HYPRE_SStructMatrixSetBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
		//	HYPRE_SStructMatrixAddToBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);

		printf("left column set\n");

		//Right column of gridpoints (Note only gets set for part1 as this portion of part0 is not on the boundary)
		bc_ilower[0] = n;
		bc_ilower[1] = 1;

		bc_iupper[0] = n;
		bc_iupper[1] = n;

		stencil_index[0] = 2;

		HYPRE_SStructMatrixSetBoxValues(A, part1, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);	
		HYPRE_SStructMatrixAddToBoxValues(A, part1,bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);

		//REMOVE THIS
		//	HYPRE_SStructMatrixSetBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
		//	HYPRE_SStructMatrixAddToBoxValues(A, part0, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);

		printf("right column set\n");

		free(valuesx);
		free(center_valuesx);
		free(valuesy);
		free(center_valuesy);

		HYPRE_SStructMatrixAssemble(A);
		HYPRE_SStructMatrixGetObject(A, (void **)&parcsr_A);

		printf("boundary conditions set\n");
	}


    /* -------------------------------------------------------------
       5. Set up SStruct Vectors for b and x
       ------------------------------------------------------------- */
	{
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
		exactsolution = calloc(nvalues* nparts, sizeof(double));

		for (j = 0; j < n; j++) {
			Y = (j + .5) * h;
			for(i = 0; i < n; i++) {
				X = (i + .5) * h;
				rhs_values[i + j * n] = -8. * PI2 * cos(2. * M_PI * X)*cos(2 * M_PI * Y);
				exactsolution[i + j * n] = 1.0 * cos(2. * M_PI * X) * cos(2. * M_PI * Y);
				x_values[i+ j * n] = 0.0;
			}

		//printf("i+j*n: %d\n", i+j*n);

		}

		//printf("%f\n", rhs_values[0]);
		HYPRE_SStructVectorSetBoxValues(b, part0, ilower, iupper, var, rhs_values);
		HYPRE_SStructVectorSetBoxValues(x, part0, ilower, iupper, var, x_values);
		HYPRE_SStructVectorSetBoxValues(x, part1, ilower, iupper, var, x_values);

		for (j = 0; j < n; j++) {
			Y = (j + .5) * h;
			for(i = 0; i < n; i++) {
				X = 1 + (i + .5) * h;
				rhs_values[i + j * n] = -8. * PI2 * cos(2. * M_PI * X)*cos(2 * M_PI * Y);
				exactsolution[i + j * n + n*n] = 1.0 * cos(2. * M_PI * X) * cos(2. * M_PI * Y);
				x_values[i+ j * n] = 0.0;
				//printf("i+j*n: %d\n", i+j*n +n*n);
			}
		}

		//printf("%f\n", rhs_values[0]);	
		HYPRE_SStructVectorSetBoxValues(b, part1, ilower, iupper, var, rhs_values);
		printf("rhs and x vectors set\n");
		free(x_values);
		free(rhs_values);

		HYPRE_SStructVectorAssemble(b);
		HYPRE_SStructVectorGetObject(b, (void **)&par_b);
		HYPRE_SStructVectorAssemble(x);
		HYPRE_SStructVectorGetObject(x, (void **)&par_x);

		HYPRE_SStructMatrixPrint("ss_halo_data/ss.initial.A", A, 0);
		HYPRE_SStructVectorPrint("ss_halo_data/ss.initial.b", b, 0);
	

    /* -------------------------------------------------------------
       6. Set up the solver and solve problem
       ------------------------------------------------------------- */
	
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
		HYPRE_ParCSRPCGSetTol(solver, 1.0e-10);
		HYPRE_ParCSRPCGSetMaxIter(solver, 1000);
		HYPRE_ParCSRPCGSetTwoNorm(solver, 1);
		HYPRE_ParCSRPCGSetLogging(solver, 3);

		// HYPRE_Solver precond;

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

		/* Writing out solution */
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
		//printf("mean = %f\n", mean);
		for (i = 1; i < totalvalues + 1; i++) {
			diff = values[i] -mean - exactsolution[i-1];
			if (isnan(diff)){
				printf("i: %d, value: %f, mean: %f, exactsolution: %f\n",i, values[i], mean, exactsolution[i-1]);		
			}
			diff2 = diff * diff;
			if (diff2 > 10){
			//printf("%+6f       %+6f     %+6f\n", values[i], exactsolution[i-1],diff);
			}
			sum += diff2;
		}	

		double L2 = sum/totalvalues;
		//printf("totalvalues = %f\n", totalvalues);
		//printf("sum = %f\n", sum);
		printf("The L2 norm is %e\n", L2);


	    /* Save the solution for GLVis visualization, see vis/glvis-ex8.sh */
	    int vis = 1;
		if (vis)
		{	
	        /* Scale indices by h */
			double T[8] = {h,0,0,h,h,0,0,h};

    	    /* Shift second part by (1,0) */
			double O[4] = {0,0,n*h,0};

			GLVis_PrintSStructGrid(grid, "vis/halo_ss.mesh", myid, T, O);
            /* Needed if solution is a ParVector */
			FILE *file;
			char filename[255];

			int nvalues = nparts*n*n;
			double *values;

            /* get the local solution */
			values = hypre_VectorData(hypre_ParVectorLocalVector(par_x));

			sprintf(filename, "%s.%06d", "vis/halo_ss.sol", myid);
			if ((file = fopen(filename, "w")) == NULL)
			{
				printf("Error: can't open output file %s\n", filename);
				MPI_Finalize();
				exit(1);
			}

            /* From vis.c */
			fprintf(file, "FiniteElementSpace\n");
            fprintf(file, "FiniteElementCollection: Local_L2_2D_P0\n");  /* cell-centered */
			fprintf(file, "VDim: 1\n");
			fprintf(file, "Ordering: 0\n\n");

            /* save solution */
			for (int i = 0; i < nvalues; i++)
			{
				fprintf(file, "%.14e\n", values[i]);
			}

			fflush(file);
			fclose(file);

			GLVis_PrintData("vis/halo_ss.data", myid, num_procs);
		}


		fflush(solution);
		fclose(solution);
	}



	HYPRE_BoomerAMGDestroy(precond);
	HYPRE_ParCSRPCGDestroy(solver);


	HYPRE_SStructMatrixDestroy(A);
	HYPRE_SStructVectorDestroy(b);
	HYPRE_SStructVectorDestroy(x);

	MPI_Finalize();

	return(0);
}




