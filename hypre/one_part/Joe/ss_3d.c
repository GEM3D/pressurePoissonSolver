/*This program solves the 3d poisson equation
 *
 * 	*insert equation*
 *
 * on a SemiStructured mesh. The code can be run in serial or in parallel.
 *
 * Compile with: make ss_3d.c
 * Execute with: mpirun -np 1 ss_3d
 */


#include <math.h>
#include <stdio.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_parcsr_ls.h"
#include "vis.c"
#include "add_to_vis.c"
#define PI2 M_PI*M_PI
#include "HYPRE_parcsr_mv.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
	//Set variables
	int i, j,k, pi, pj, pk;
	int myid, num_procs;
	int Nx, Ny, Nz, nx, ny, nz;
	int entry, nvalues, nvaluesx, nvaluesy, nvaluesz;
	int solver_id, vis;
	int num_iterations;
	int ilower[3], iupper[3],stencil_indices[7];
	int bc_ilower[3], bc_iupper[3], center_index[1], stencil_index[1];
	int part = 0;
	int nparts = 1;
	int nvars = 1;
	int var = 0;
	int nentries = 7;
	int index = 1;
	
	int object_type = HYPRE_PARCSR;


	double hx, hy, hz;
	double X, Y, Z;
	double hx2inv, hy2inv, hz2inv;
	double *values;
	double *valuesx, *valuesy, *valuesz, *center_valuesx, *center_valuesy, *center_valuesz;
	double *rhs_values, *x_values, *exactsolution;
	double final_res_norm, t1, t2;

	HYPRE_SStructGraph graph;
	HYPRE_SStructGrid grid;
	HYPRE_SStructStencil stencil;
	HYPRE_SStructMatrix A;
	HYPRE_ParCSRMatrix par_A;
	HYPRE_SStructVector b;
	HYPRE_SStructVector x;
	HYPRE_ParVector par_b;
	HYPRE_ParVector par_x;
printf("variables set\n");	
	//Initialize MPi
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
printf("mpi setup\n");
	//Set default problem parameters
	nx = 100;
	ny = 100;
	nz = 100;
	Nx = 1;
	Ny = 1;
	Nz = 1;
	solver_id = 0;
	vis = 1;
printf("default parameters set\n");
	nvalues = nentries * nx * ny;
	hx = 1.0 / (Nx * nx);
	hy = 1.0 / (Ny * ny);
	hz = 1.0 / (Nz * nz);
printf("gird spacing calculated\n");
	//Set the processors in a grid
	pi = 0;
	pj = 0;
	pk = 0;

	//Determine each processors part of the grid
	ilower[0] = 0;
	ilower[1] = 0;
	ilower[2] = 0;

	iupper[0] = nx + ilower[0] - 1;
	iupper[1] = ny + ilower[1] - 1;
	iupper[2] = nz + ilower[2] - 1;
printf("processor grid and division set\n");
	//Create and setup a 3d grid
	HYPRE_SStructGridCreate(MPI_COMM_WORLD, 3, nparts, &grid);
	HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);
printf("3d grid setup\n");
	//Set variable types to cell centered
	HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};

	for(i = 0; i < nparts; i++){
		HYPRE_SStructGridSetVariables(grid, i, nvars, vartypes);
	}
printf("variable types set\n");
	//Assemble the grid
	HYPRE_SStructGridAssemble(grid);
	
	//Set up the 3d, 7pt stencil
	HYPRE_SStructStencilCreate(3, 7, &stencil);
	int offsets[7][3] = {{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
	
	for (entry = 0; entry < 5; entry++) {
		HYPRE_SStructStencilSetEntry(stencil,entry, offsets[entry], var);	
	}
printf("stencil created\n");
	//Create the graph
	HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);
	
	//Change the graph's type to be sparse
	HYPRE_SStructGraphSetObjectType(graph, object_type);
	
	//Apply the 3d 7pt stencil to the graph
	HYPRE_SStructGraphSetStencil(graph, part, var, stencil);
	
	//Assemble the graph
	HYPRE_SStructGraphAssemble(graph);

	//Setup the matrix
	HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
	HYPRE_SStructMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_SStructMatrixInitialize(A);
	
	//Set the matrix coefficients
	values = calloc(nvalues, sizeof(double));
	for (i = 0; i < nentries; i++){
		stencil_indices[i] = i;
	}

	//Account for differences in hx, hy, and hz;
	hx2inv = 1 / (hx * hx);
	hy2inv = 1 / (hy * hy);
	hz2inv = 1 / (hz * hz);

	for(i = 0; i < nvalues; i += nentries) {
		values[i] = -2.0 * hx2inv -2.0 * hy2inv;
		values[i+1] = hx2inv;
		values[i+2] = hx2inv;
		values[i+3] = hy2inv;
		values[i+4] = hy2inv;
		values[i+5] = hz2inv;
		values[i+6] = hz2inv;
	}

	//Set the stencil values
	HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, var, nentries, stencil_indices, values);
printf("Stencils set\n");
	//Set the Neumann boundary condition
	index = 1;
	nvaluesx = index * nx * ny;
	nvaluesy = index * nx * nz;
	nvaluesz = index * nz *ny;
	center_index[0] = 0;

	valuesx = calloc(nvaluesx, sizeof(double));
	valuesy = calloc(nvaluesy, sizeof(double));
	valuesz = calloc(nvaluesz, sizeof(double));

	center_valuesx = calloc(nvaluesx, sizeof(double));
	center_valuesy = calloc(nvaluesy, sizeof(double));
	center_valuesz = calloc(nvaluesz, sizeof(double));

	for(i = 0; i < nvaluesx; i++){
		valuesx[i] = 0.0;
		center_valuesx[i] = hx2inv;
	}

	for(i = 0; i < nvaluesy; i++){
		valuesy[i] = 0.0;
		center_valuesy[i] = hy2inv;
	}

	for(i = 0; i < nvaluesz; i++){
		valuesz[i] = 0.0;
		center_valuesz[i] = hz2inv;
	}

printf("bc memory allocated\n");
	//Recall: pi, pj, and pk describe the position in the precessor grid
	//Bottom face
	if(pj == 0){
		bc_ilower[0] = pi * nx;
		bc_ilower[1] = pj * ny;
		bc_ilower[2] = pk * nz;
	
		bc_iupper[0] = bc_ilower[0] + nx -1;
		bc_iupper[1] = bc_ilower[1];
		bc_iupper[2] = bc_ilower[2] + nz -1;
		//Set the index corresponding to the bottom of the stencil
		stencil_index[0] = 3;

		//Sets one boundary condition to be Dirichlet (allowing an exact solution to be obtained)
		if (myid == 0){
			center_valuesy[1] = -hy2inv;
		}

printf("bc_lower: [%i,%i,%i]\n", bc_ilower[0], bc_ilower[1], bc_ilower[2]);
printf("bc_upper: [%i,%i,%i]\n", bc_iupper[0], bc_iupper[1], bc_iupper[2]);
printf("bc bot start\n");
		//Sets the boundary condition on the bottom face of the domain. Set outsides to 0 and hx2inv to the center entry
		HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesy);
printf("bc bot box vals set\n");
		HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesy);
printf("bc bot box add\n");
		//Changes back to the full Neumann condition
		center_valuesy[1] = hy2inv;

	}	
printf("bottom face set\n");	
	//Top face
	if (pj == Ny -1){
		bc_ilower[0] = pi * nx;
		bc_ilower[1] = pj * ny + ny - 1;
		bc_ilower[2] = pk * nz;

		bc_iupper[0] = bc_ilower[0] + nx - 1;
		bc_iupper[1] = bc_ilower[1];
		bc_iupper[2] = bc_ilower[2] + nz -1;

		//Set the index corresponding to the top of the stencil
		stencil_index[0] = 4;

		//Enforces the Neumann condition on the top face of the domain
		HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesy);
		HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesy);
	}
printf("top face set\n");
	//Left face
	if (pi == 0){
		bc_ilower[0] = pi * nx;
		bc_ilower[1] = pj * ny;
		bc_ilower[2] = pk * nz;

		bc_iupper[0] = bc_ilower[0];
		bc_iupper[1] = bc_ilower[1] + ny -1;
		bc_iupper[2] = bc_ilower[2] + nz -1;

		//Set the index corresponding to the left of the stencil
		stencil_index[0] = 1;

		//Enforces the Neumann condition on the left face of the domain
		HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
		HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);
	}
printf("left face set\n");
	//Right face
	if (pi == Nx -1){
		bc_ilower[0] = pi * nx + nx -1;
		bc_ilower[1] = pj * ny;
		bc_ilower[2] = pk * nz;

		bc_iupper[0] = bc_ilower[0];
		bc_iupper[1] = bc_ilower[1] + ny -1;
		bc_iupper[2] = bc_ilower[2] + nz -1;

		//Set the index corresponding to the right of the stencil
		stencil_index[0] = 2;

		//Enforces the Neumann condition on the right face of the domain
		HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesx);
		HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);
	}
printf("right face set\n");
	//Back face
	if (pk == Nz-1){
		bc_ilower[0] = pi * nx;
		bc_ilower[1] = pj * ny; 
		bc_ilower[2] = pk * nz + nz -1;

		bc_iupper[0] = bc_ilower[0] + nx -1;
		bc_iupper[1] = bc_ilower[1] + ny -1;
		bc_iupper[2] = bc_ilower[2];

		//Set the index corresponding to the back of the stencil
		stencil_index[0] = 5;

		//Enforces the Neumann condition on the back face of the domain
		HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesz);
		HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesz);
	}
printf("back face set\n");
	//Front face
        if (pk == 0){
                bc_ilower[0] = pi * nx;
                bc_ilower[1] = pj * ny;
                bc_ilower[2] = pk * nz;

                bc_iupper[0] = bc_ilower[0] + nx -1;
                bc_iupper[1] = bc_ilower[1] + ny -1;
                bc_iupper[2] = bc_ilower[2];

		//Set the index correspoinding to the front of the stencil
		stencil_index[0] = 6;

		//Enforces the Neumann condition of the front face of the domain
		HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, stencil_index, valuesz);
		HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesz);
	}
printf("front face set\n");

	//Free memory used to enforce boundary conditions
	free(valuesx);
	free(valuesy);
	free(valuesz);
	free(center_valuesx);
	free(center_valuesy);
	free(center_valuesz);
printf("face boundary conditions set\n");
	//Assemble the Matrix
	HYPRE_SStructMatrixAssemble(A);

	//Fills par_A with the sparse version of A
	HYPRE_SStructMatrixGetObject(A, (void **)&par_A);
	
	//Create the right hand side and solution vectors
	HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &b);
	HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &x);

	//Set the right hand side and solution vectors to be sparse
	HYPRE_SStructVectorSetObjectType(b, HYPRE_PARCSR);
	HYPRE_SStructVectorSetObjectType(x, HYPRE_PARCSR);
	

	//Initialize the right hand side and solution vectors
	HYPRE_SStructVectorInitialize(b);
	HYPRE_SStructVectorInitialize(x);

	//Allocate memory for the right hand side, exact solution, and initial guess values
        nvalues = nx * ny * nz;
	
	rhs_values = calloc(nvalues, sizeof(double));
	x_values = calloc(nvalues, sizeof(double));
	exactsolution = calloc(nvalues, sizeof(double));
	
printf("begin initial value calculation\n");
	//Calculate the rhs and exact solution at each point
	for(k = 0; k < nz; k++){
		Z = (ilower[2] + k + .5) * hz;
		for(j = 0; j < ny; j++){
			Y = (ilower[1] + j + .5) * hy;
			for(i = 0; j < nx; i++){
				X = (ilower[0] + i + .5) * hx;
				//Change this to a real calculation
				rhs_values[i+j*nx+k*nx*ny] = (X+Y+Z)*0. + 1.;
				exactsolution[i+j*nx+k*nx*ny] = 1.;
				x_values[i+j*nx+k*nx*ny] = 0.0;
				
			}
		}
	}

printf("initial values calculated\n");


	//Set values for the rhs and exact solution vectors
	HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, var, rhs_values);
	HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, var, x_values);

	//Free memory used in calculation of right hand side and x vectors
	free(rhs_values);
	free(x_values);		

	//Assemble right hand side and x vectors
	HYPRE_SStructVectorAssemble(b);
	HYPRE_SStructVectorAssemble(x);

	//Place the sparse version of b and x into par_b and par_x respectively
	HYPRE_SStructVectorGetObject(b, (void **)&par_b);
	HYPRE_SStructVectorGetObject(x, (void **)&par_x);

	//Print out the initial Matrix, A and right hand side vector b
	HYPRE_SStructMatrixPrint("3d_ss_data/ss.initial.A", A, 0);
	HYPRE_SStructVectorPrint("3d_ss_data/ss.initial.b", b, 0);

	//Solve the system
	//Boomer AMG
	if (solver_id == 0) {
		//Create GMRES solver
		HYPRE_Solver gmres_solver;

		//Set GMRES solver parameters
		HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &gmres_solver);
		HYPRE_GMRESSetMaxIter(gmres_solver, 1000);
		HYPRE_GMRESSetTol(gmres_solver, 1e-12);
		HYPRE_GMRESSetPrintLevel(gmres_solver, 0);
		HYPRE_GMRESSetLogging(gmres_solver, 0);

		//Create AMG preconditioner
		HYPRE_Solver precond;

		//Set AMG preconditioner parameters
		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
        	HYPRE_BoomerAMGSetCoarsenType(precond, 6);
        	HYPRE_BoomerAMGSetTol(precond, 1.0e-3);
        	HYPRE_BoomerAMGSetPrintLevel(precond, 0);
        	HYPRE_BoomerAMGSetMaxIter(precond, 1);
        	HYPRE_BoomerAMGSetNumSweeps(precond, 2);
		HYPRE_BoomerAMGSetRelaxType(precond, 0);
        	//HYPRE_BoomerAMGSetRelaxWeight(precond, .3);
                HYPRE_BoomerAMGSetRelaxWt(precond, .3);
                
		//Setup and solve the system
		HYPRE_ParCSRGMRESSetPrecond(gmres_solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
        	HYPRE_ParCSRGMRESSetup(gmres_solver, par_A, par_b, par_x);
        	t1 = MPI_Wtime();
        	HYPRE_ParCSRGMRESSolve(gmres_solver, par_A, par_b, par_x);
        	t2 = MPI_Wtime();	

		//Get run info
		HYPRE_ParCSRGMRESGetNumIterations(gmres_solver, &num_iterations);
	        HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(gmres_solver, &final_res_norm);

		//Cleanup
	        HYPRE_ParCSRGMRESDestroy(gmres_solver);
	        HYPRE_BoomerAMGDestroy(precond);

	}

	MPI_Barrier(MPI_COMM_WORLD);


	//Uses MPI to find the maximum time taken to solve the equation
	double *time, *timemax;
	time = calloc(1, sizeof(double));
  	time[0] = t2-t1;
  	timemax = calloc(1, sizeof(double));
  	int root = 0;
  	int count = 1;
  	MPI_Reduce(time, timemax, count, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

	//Prints the maximum time taken to solve the problem
  	if (myid == 0) {
        	printf("\n");
        	printf("Elapsed Wall Time %f\n", timemax[0]);
  	}


	//Saves the solution for visualization and L2 norm calculation
	if (vis) {
		//Creates a file for the solution
		FILE *file;
    		FILE *solution;

		//Names the solution file
    		char filename[255];
    		char solutionfile[255];

    		sprintf(solutionfile, "%s.%06d.%d", "3d.ss.final.x", myid, myid);

		//Determines the number of values the solution file needs to hold
    		int nvalues = nx * ny * nz;

		int root = 0;
		int count = 1;
    		double sum = 0;
		double mean = 0;
		double totalvalues = nvalues * num_procs * 1.0;
    		double diff, diff2, L2;
		double *sendbuffer, *recvbuffer;
         	double *values = calloc(nvalues+1, sizeof(double));
    
		//Opens the solution file, prints an error message if the file cannot be opened
		if ((solution = fopen(solutionfile, "rt")) == NULL) {
		  	printf("Error: can't open output file %s\n", solutionfile);
      			MPI_Finalize();
      			exit(1);
    		}

		//Reads the solution file into the array, values and calculates the mean solution value
		for (i = 0; i < nvalues + 1; i++) {
		        fscanf(solution, "%lf", &values[i]);
      			if (i != 0){
      				mean += values[i];
      			}	
    		}	

		//Closes the solution file
		fflush(solution);
		fclose(solution);

    		MPI_Barrier(MPI_COMM_WORLD);
		mean /= nvalues;

		//Calculates the sum of the squared differences between the exact and numerical solutions on each processor
		for (i = 1; i < nvalues + 1; i++) {
		        diff = values[i] -mean - exactsolution[i-1];
     			diff2 = diff * diff;
      			sum += diff2;
      		}			

		//Makes sure all processors are caught up
		MPI_Barrier(MPI_COMM_WORLD);

		//Passes the individual sums to processor 0 and sums them
		sendbuffer = calloc(1,  sizeof(double));
    		sendbuffer[0] = sum;
    		recvbuffer = calloc(1, sizeof(double));
    		MPI_Reduce(sendbuffer, recvbuffer, count, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    		MPI_Barrier(MPI_COMM_WORLD);

		//Calculates the L2 norm
		if (myid == 0){
			L2 = sqrt(recvbuffer[0] / totalvalues);
			printf("The L2 norm is %e\n", L2);
		}

		//Names the visualization file
		sprintf(filename, "%s.%06d", "vis/3d.ss.sol", myid);
		
		//Opens the visualization file, and if the file cannot be open an error is returned
		if ((file = fopen(filename, "w")) == NULL) {
      			printf("Error: can't open output file %s\n", filename);
      			MPI_Finalize();
      			exit(1);
    		}

		//Print the values to the visualization file
		for (i = 1; i < nvalues+1; i++) {
			fprintf(file, "%.14e\n", values[i]);
		}

		//Vis Cleanup
		fflush(file);
		fclose(file);

		free(values);
		free(exactsolution);
		free(sendbuffer);
		free(recvbuffer);

		//Save global finite element mesh
		if (myid == 0){
			GLVis_PrintGlobalMesh("vis/ss.mesh", Nx, Ny, nx, ny, hx, hy);
		}
	}

	//Clean up
	HYPRE_SStructMatrixDestroy(A);
	HYPRE_SStructVectorDestroy(b);
	HYPRE_SStructVectorDestroy(x);

	//Finalize MPI
	MPI_Finalize();

	return(0);
}
	










