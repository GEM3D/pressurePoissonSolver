/* This program uses the SStruct format to set up a poisson equation, converts
 * the matrices and vectors to a sparse format and solves using AMGBOOMER
 *
 * The problem to be solved is the poisson equation where u =
 *sin(pi*x)*sin(pi*y)
 * f = -2pi^2*sin(2*pi*x)*(2*pi*y) with dirchlet boundary conditions
 *
 * Compile with: make ExactIJNeumann
 * Execute with: mpirun -np 4 ExactIJNeumann -nx 100 -ny 100 -NX 2 -solver 0 -vis
 *
 * Create by: Joseph McNeal
 */

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
  // printf("start\n");
  int i, j, pi, pj;
  int myid, num_procs;
  int Nx, Ny, nx, ny;
  //nx and ny are number of mesh cells
  int solver_id;
  int vis;
  int ilower[2], iupper[2];
  int nparts = 1;
  int part = 0;

  double hx, hy;

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
  // printf("variables initialized\n");
   
  //Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  // printf("mpi initialized\n");
   
  // Set default problem parameters
  nx = 100;
  ny = 100;
  Nx = 1;
  solver_id = 0;
  vis = 0;
  // printf("defaults set\n");
   
  // Parse command line

  int arg_index = 0;

  while (arg_index < argc) {
    if (strcmp(argv[arg_index], "-nx") == 0) {
      arg_index++;
      nx = atoi(argv[arg_index++]);
      // printf("nx set\n");
    } else if (strcmp(argv[arg_index], "-ny") == 0) {
      arg_index++;
      ny = atoi(argv[arg_index++]);
      // printf("ny set\n");
    } else if (strcmp(argv[arg_index], "-Nx") == 0) {
      arg_index++;
      Nx = atoi(argv[arg_index++]);
      // printf("Nx set\n");
    } else if (strcmp(argv[arg_index], "-solver") == 0) {
      arg_index++;
      solver_id = atoi(argv[arg_index++]);
      // printf("solver set\n");
    } else if (strcmp(argv[arg_index], "-vis") == 0) {
      arg_index++;
      vis = 1;
      // printf("print solution\n");
    } else {
      arg_index++;
      // printf("plus plus\n");
    }
  }

  // printf("command line parsed\n");
   
  // Set up processors in grid 
  Ny = num_procs / Nx;
  hx = 1.0 / (Nx * nx);
  hy = 1.0 / (Ny * ny);

  if (Ny <= Nx) {
    pj = myid / Nx;
    pi = myid - pj * Ny;
  } else {
    pi = myid / Ny;
    pj = myid - pi * Nx;
  }
  // printf("processors assigned\n");

  // Determine each processor's piece of the grid 
  ilower[0] = pi * nx;
  ilower[1] = pj * ny;

  iupper[0] = ilower[0] + nx - 1;
  iupper[1] = ilower[1] + ny - 1;

  // printf("upper and lowers assigned\n");
  HYPRE_SStructGridCreate(MPI_COMM_WORLD, 2, nparts, &grid);
  // printf("grid created\n");
  HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);
  // printf("grid set\n");
  
  int nvars = 1;

  HYPRE_SStructVariable vartypes[1] = { HYPRE_SSTRUCT_VARIABLE_CELL };
  // printf("variable type set\n");
  for (i = 0; i < nparts; i++) {
    HYPRE_SStructGridSetVariables(grid, i, nvars, vartypes);
  }
  // printf("variable type set\n");
  HYPRE_SStructGridAssemble(grid);
  // printf("grid assembled\n");
  
  // Set up the stencil 
  HYPRE_SStructStencilCreate(2, 5, &stencil);
  int entry;
  int offsets[5][2] = { { 0, 0 }, { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } };
  int var = 0;

  for (entry = 0; entry < 5; entry++) {
    HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], var);
  }

  // Create the graph 
  HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);

  int object_type = HYPRE_PARCSR;

  HYPRE_SStructGraphSetObjectType(graph, object_type);

  HYPRE_SStructGraphSetStencil(graph, part, var, stencil);
  HYPRE_SStructGraphAssemble(graph);

  int nentries = 5;
  int nvalues = nentries * nx * ny;
  double *values;
  int stencil_indices[5];

  // Setup the matrix 
  HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
  HYPRE_SStructMatrixSetObjectType(A, HYPRE_PARCSR);
  HYPRE_SStructMatrixInitialize(A);

  // Set the matrix coefficients 
  values = calloc(nvalues, sizeof(double));
  for (j = 0; j < nentries; j++) {
    stencil_indices[j] = j;
  }

  // Account for differences in hx and hy 
  double hx2 = hx * hx;
  double hy2 = hy * hy;
  double hx2inv = 1 / hx2;
  double hy2inv = 1 / hy2;

  for (i = 0; i < nvalues; i += nentries) {
    values[i] = -2.0 * hx2inv - 2.0 * hy2inv;
    values[i + 1] = hx2inv;
    values[i + 2] = hx2inv;
    values[i + 3] = hy2inv;
    values[i + 4] = hy2inv;
  }

  HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, var, nentries,
                                  stencil_indices, values);

  // Set the Neumann boundary condition
  {
    int bc_ilower[2];
    int bc_iupper[2];
    int nentries = 1;
    int nvaluesx = nentries * nx; /*  number of stencil entries times the length
                                      of one side of my grid box */
    int nvaluesy = nentries * ny;
    double *valuesx, *valuesy, *center_valuesx, *center_valuesy;
    int stencil_indices[1];
    int center_index[1];
    center_index[0] = 0;

    valuesx = calloc(nvaluesx, sizeof(double));
    center_valuesx = calloc(nvaluesx, sizeof(double));
    for (j = 0; j < nvaluesx; j++){
      valuesx[j] = 0.0;
      center_valuesx[j] = hx2inv;
    }

    valuesy = calloc(nvaluesy, sizeof(double));
    center_valuesy = calloc(nvaluesy, sizeof(double));
    for (j = 0; j < nvaluesy; j++){
      valuesy[j] = 0.0;
      center_valuesy[j] = hy2inv;
     }
 

    // Recall: pi and pj describe position in the processor grid 
    if (pj == 0) {
      // Bottom row of grid points 
      bc_ilower[0] = pi * nx;
      bc_ilower[1] = pj * ny;

      bc_iupper[0] = bc_ilower[0] + nx - 1;
      bc_iupper[1] = bc_ilower[1];

      stencil_indices[0] = 3;

      //Sets one of the boundary conditions to the dirichlet condition
      if (myid == 0){
	center_valuesy[1] = -hy2inv;
      }
      
      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var,
                                      nentries, stencil_indices, valuesy);
      HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesy);
      center_valuesy[1] = hy2inv;
      // printf("Bottom row zeroed\n");
    }

    if (pj == Ny - 1) {
      // upper row of grid points 
      bc_ilower[0] = pi * nx;
      bc_ilower[1] = pj * ny + ny - 1;

      bc_iupper[0] = bc_ilower[0] + ny - 1;
      bc_iupper[1] = bc_ilower[1];

      stencil_indices[0] = 4;

      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var,
                                      nentries, stencil_indices, valuesy);
      HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesy);
      // printf("Upper row zeroed\n");
    }
    if (pi == 0) {
      // Left row of grid points 
      bc_ilower[0] = pi * nx;
      bc_ilower[1] = pj * ny;

      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + ny - 1;

      stencil_indices[0] = 1;

      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var,
                                      nentries, stencil_indices, valuesx);
      HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx); 
      // printf("left collumn zeroed\n");
    }

    if (pi == Nx - 1) {
      // Right row of grid points 
      bc_ilower[0] = pi * nx + nx - 1;
      bc_ilower[1] = pj * ny;

      bc_iupper[0] = bc_ilower[0];
      bc_iupper[1] = bc_ilower[1] + ny - 1;

      stencil_indices[0] = 2;

      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var,
                                      nentries, stencil_indices, valuesx);
      HYPRE_SStructMatrixAddToBoxValues(A, part, bc_ilower, bc_iupper, var, nentries, center_index, center_valuesx);
      // printf("right collumn zeroed\n");
    }
    free(valuesx);
    free(center_valuesx);
    free(valuesy);
    free(center_valuesy);
  }

  HYPRE_SStructMatrixAssemble(A);
  HYPRE_SStructMatrixGetObject(A, (void **)&parcsr_A);

  // Create the right hand side and solution vectors 
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &b);
  HYPRE_SStructVectorSetObjectType(b, HYPRE_PARCSR);
  HYPRE_SStructVectorInitialize(b);

  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &x);
  HYPRE_SStructVectorSetObjectType(x, HYPRE_PARCSR);
  HYPRE_SStructVectorInitialize(x);

  // Set the right hand side, exact solution, and initial guess values 
  double *rhs_values, *x_values, *exactsolution;
  nvalues = nx * ny;

  rhs_values = calloc(nvalues, sizeof(double));
  x_values = calloc(nvalues, sizeof(double));
  exactsolution = calloc(nvalues, sizeof(double));
  // printf("RHS calculation begin\n");
    for (j = 0; j < ny; j++) {
         double y =  (ilower[1] + j + .5) * hy;
    for (i = 0; i < nx; i++) {
         double x = (ilower[0] + i + .5) * hx;
      rhs_values[i + j * nx] =
          -8. * PI2 * cos(2. * M_PI * x) * cos(2. * M_PI * y);
      exactsolution[i + j * nx] = 1.0 * cos(2. * M_PI * x) * cos(2. * M_PI * y);
      x_values[i + j * nx] = 0.0;
    }
  }
  // printf("RHS calculation ended\n");
  HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, var, rhs_values);
  HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, var, x_values);


  //printf("%d  %f\n", (ny*nx-1), exactsolution[ny*nx-1]);
  

  free(x_values);
  free(rhs_values);

  HYPRE_SStructVectorAssemble(b);
  // HYPRE_SStructVectorPrint("ss.initial.b", b, 0);
  HYPRE_SStructVectorGetObject(b, (void **)&par_b);
  HYPRE_SStructVectorAssemble(x);
  HYPRE_SStructVectorGetObject(x, (void **)&par_x);

  // Print out initial vectors 
  HYPRE_SStructMatrixPrint("SStructExact_data/ss.initial.A", A, 0);
  HYPRE_SStructVectorPrint("SStructExact_data/ss.initial.b", b, 0);
  // printf("About to solve\n");

  int num_iterations;
  double final_res_norm, t1, t2;


  // Select a solver 
  // BoomerAMG
  if (solver_id == 0) {

    // Create AMG solver 
    HYPRE_BoomerAMGCreate(&solver);

    // Set solver parameters 
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);
    HYPRE_BoomerAMGSetCoarsenType(solver, 6);
    HYPRE_BoomerAMGSetRelaxType(solver, 3);
    HYPRE_BoomerAMGSetNumSweeps(solver, 1);
    HYPRE_BoomerAMGSetMaxLevels(solver, 20);
    HYPRE_BoomerAMGSetTol(solver, 1.0e-12);

    //HYPRE_BoomerAMGSetMaxIter(solver, 100);
    // Set up and solve 
    HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
    t1 = MPI_Wtime();
    HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);
    t2 = MPI_Wtime();
    //printf("solution obtained\n");
    
    
    //Destroy solver
    HYPRE_BoomerAMGDestroy(solver);
  
    // Run information 
    HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
  }  
  
  //BoomerAMG preconditioned PCG
  if (solver_id == 1) {
	//Create PCG solver
	HYPRE_Solver pcgsolver;

	//Set PCG solver parameters
	HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &pcgsolver);
	HYPRE_ParCSRPCGSetTol(pcgsolver, 1.0e-12);
	HYPRE_ParCSRPCGSetMaxIter(pcgsolver, 500);
	HYPRE_ParCSRPCGSetTwoNorm(pcgsolver, 1);
	HYPRE_ParCSRPCGSetLogging(pcgsolver, 3);

	//Create AMG preconditioner
	HYPRE_Solver precond;
	
	//Set AMG preconditioner parameters
	HYPRE_BoomerAMGCreate(&precond);
	HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
	HYPRE_BoomerAMGSetCoarsenType(precond, 6);
	HYPRE_BoomerAMGSetTol(precond, 1.0e-12);
	HYPRE_BoomerAMGSetPrintLevel(precond, 1);
        HYPRE_BoomerAMGSetPrintLevel(precond, 0);
        HYPRE_BoomerAMGSetMaxIter(precond, 5);
        HYPRE_BoomerAMGSetNumSweeps(precond, 2);
        HYPRE_BoomerAMGSetRelaxType(precond, 0);
        //HYPRE_BoomerAMGSetRelaxWeight(precond, .3);
        HYPRE_BoomerAMGSetRelaxWt(precond, .3);
                

	
  	//Setup and solve the system
  	HYPRE_ParCSRPCGSetPrecond(pcgsolver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
	HYPRE_ParCSRPCGSetup(pcgsolver, parcsr_A, par_b, par_x);
	t1 = MPI_Wtime();
	HYPRE_ParCSRPCGSolve(pcgsolver, parcsr_A, par_b, par_x);
	t2 = MPI_Wtime();

	//Get run info
	HYPRE_ParCSRPCGGetNumIterations(pcgsolver, &num_iterations);
	HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(pcgsolver, &final_res_norm);

	//Cleanup
	HYPRE_BoomerAMGDestroy(precond);
	HYPRE_ParCSRPCGDestroy(pcgsolver);
  }

  //Boomer AMG Preconditioned GMRES
  if (solver_id == 2){
	//Create GMRES solver
	HYPRE_Solver gmres_solver;

	//Set GMRES solver parameters
	HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &gmres_solver);
	HYPRE_GMRESSetMaxIter(gmres_solver, 1000);
	//HYPRE_GMRESSetAbsoluteTol(gmres_solver, 1.0e-12);
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
	HYPRE_ParCSRGMRESSetup(gmres_solver, parcsr_A, par_b, par_x);
	t1 = MPI_Wtime();
	HYPRE_ParCSRGMRESSolve(gmres_solver, parcsr_A, par_b, par_x);
	t2 = MPI_Wtime();

	//Get run info
	HYPRE_ParCSRGMRESGetNumIterations(gmres_solver, &num_iterations);
	HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(gmres_solver, &final_res_norm);

	//Cleanup
	HYPRE_ParCSRGMRESDestroy(gmres_solver);
	HYPRE_BoomerAMGDestroy(precond);
  }

  //Euclid (ILU) preconditioned GMRES
  if (solver_id == 3){
	//Create GMRES solver
	HYPRE_Solver gmres_solver;
	HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &gmres_solver);

	//Set the GMRES solver parameters
	HYPRE_GMRESSetMaxIter(gmres_solver,100);
	HYPRE_GMRESSetTol(gmres_solver, 1e-12);
	HYPRE_GMRESSetPrintLevel(gmres_solver, 3);
	HYPRE_GMRESSetLogging(gmres_solver, 0);	
	HYPRE_GMRESSetKDim(gmres_solver, 1);

	//Create Euclid precondtioner
	HYPRE_Solver euclid_precond;
	HYPRE_EuclidCreate(MPI_COMM_WORLD, &euclid_precond);
	//HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &euclid_precond);

	//Set the euclid precondioner parameters
	HYPRE_EuclidSetLevel(euclid_precond, 1);
	//HYPRE_EuclidSetBJ(euclid_precond, 0);
	HYPRE_ParaSailsSetParams(euclid_precond, .01, 2);
	//HYPRE_ParaSailsSetSym(euclid_precond, 2);


	//Setup and solve the system
	HYPRE_ParCSRGMRESSetPrecond(gmres_solver, HYPRE_EuclidSolve, HYPRE_EuclidSetup, euclid_precond);
	//HYPRE_ParCSRGMRESSetPrecond(gmres_solver, HYPRE_ParaSailsSolve, HYPRE_ParaSailsSetup, euclid_precond);
	HYPRE_ParCSRGMRESSetup(gmres_solver, parcsr_A, par_b, par_x);
	t1 = MPI_Wtime();
	HYPRE_ParCSRGMRESSolve(gmres_solver, parcsr_A, par_b, par_x);
	t2 = MPI_Wtime();

	//Get run info
	HYPRE_ParCSRGMRESGetNumIterations(gmres_solver, &num_iterations);
	HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(gmres_solver, &final_res_norm);

	//Cleanup
	HYPRE_ParCSRGMRESDestroy(gmres_solver);
	HYPRE_EuclidDestroy(euclid_precond);
  }



  //Doesn't work because the ParCSR solvers want Par CSR preconditioners
  //SMG (GMG) preconditioned PCG
  /*
  if (solver_id == 4){
	//Set PCG parameters
	HYPRE_Solver pcg_solver;
	HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &pcg_solver);
        HYPRE_ParCSRPCGSetTol(pcg_solver, 1.0e-12);
        HYPRE_ParCSRPCGSetMaxIter(pcg_solver, 500);
        HYPRE_ParCSRPCGSetTwoNorm(pcg_solver, 1);
        HYPRE_ParCSRPCGSetLogging(pcg_solver, 3);
	
	//Set SMG parameters
        HYPRE_Solver smg_precond;
	HYPRE_StructSMGCreate(MPI_COMM_WORLD, &smg_precond);
	HYPRE_StructSMGSetMaxIter(smg_precond, 1);
	HYPRE_StructSMGSetTol(smg_precond, 1.0e-12);	
	//Setup and solve the system
	HYPRE_ParCSRPCGSetPrecond(pcg_solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, smg_precond);
        HYPRE_ParCSRPCGSetup(pcg_solver, parcsr_A, par_b, par_x);
        t1 = MPI_Wtime();
        HYPRE_ParCSRPCGSolve(pcg_solver, parcsr_A, par_b, par_x);
        t2 = MPI_Wtime();

	//Get run info
	HYPRE_ParCSRPCGGetNumIterations(pcg_solver, &num_iterations);
	HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(pcg_solver, &final_res_norm);
	
	//cleanup
	HYPRE_ParCSRPCGDestroy(pcg_solver);
	}
     */

  //AMG preconditioned BiCGSTAB
  if (solver_id == 5){

	//Set BiCGSTAB parameters
	HYPRE_Solver bicg_solver;
	HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &bicg_solver);
	HYPRE_ParCSRBiCGSTABSetTol(bicg_solver, 1.0e-12);
	HYPRE_ParCSRBiCGSTABSetPrintLevel(bicg_solver, 3);
	HYPRE_ParCSRBiCGSTABSetLogging(bicg_solver, 0);
	HYPRE_ParCSRBiCGSTABSetMaxIter(bicg_solver, 100);

	//Set AMG parameters
	HYPRE_Solver amg_precond;
	HYPRE_BoomerAMGCreate(&amg_precond);
        HYPRE_BoomerAMGSetStrongThreshold(amg_precond, .25);
        HYPRE_BoomerAMGSetCoarsenType(amg_precond, 6);
        HYPRE_BoomerAMGSetTol(amg_precond, 1.0e-12);
        HYPRE_BoomerAMGSetPrintLevel(amg_precond, 1);
        HYPRE_BoomerAMGSetPrintLevel(amg_precond, 0);
        HYPRE_BoomerAMGSetMaxIter(amg_precond, 5);
        HYPRE_BoomerAMGSetNumSweeps(amg_precond, 2);
        HYPRE_BoomerAMGSetRelaxType(amg_precond, 0);
        //HYPRE_BoomerAMGSetRelaxWeight(amg_precond, .3);
        HYPRE_BoomerAMGSetRelaxWt(amg_precond, .3);
	
	//Setup and solve the system
	HYPRE_ParCSRBiCGSTABSetPrecond(bicg_solver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, amg_precond);
	HYPRE_ParCSRBiCGSTABSetup(bicg_solver, parcsr_A, par_b, par_x);
	t1 = MPI_Wtime();
	HYPRE_ParCSRBiCGSTABSolve(bicg_solver, parcsr_A, par_b, par_x);
	t2 = MPI_Wtime();
	//Get run info
	HYPRE_ParCSRBiCGSTABGetNumIterations(bicg_solver, &num_iterations);
	HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(bicg_solver, &final_res_norm);
	//Cleanup
	HYPRE_ParCSRBiCGSTABDestroy(bicg_solver);
	HYPRE_BoomerAMGDestroy(amg_precond);
	}


  //Euclid precondtioned BiCGSTAB
  if (solver_id == 6){
	}


    if (myid == 0) {
      printf("\n");
      printf("Iterations = %d\n", num_iterations);
      printf("Final Relative Residual Norm = %e\n", final_res_norm);
      printf("\n");
    }

    // HYPRE_ParVectorPrint(par_x, "SStructExact_data/ss.final.x");
    char solutionfile[255];
    sprintf(solutionfile, "%s.%06d", "ss.final.x", myid);
    HYPRE_ParVectorPrint(par_x, solutionfile);
    
  MPI_Barrier(MPI_COMM_WORLD);

  double *time, *timemax;
  time = calloc(1, sizeof(double));
  time[0] = t2-t1;
  timemax = calloc(1, sizeof(double));
  int root = 0;
  int count = 1;
  MPI_Reduce(time, timemax, count, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

  if (myid == 0) {
	printf("\n");
	printf("Elapsed Wall Time %f\n", timemax[0]);
  }


  // Save the solution for visualization and L2 norm calculation  
  if (vis) {
    FILE *file;
    FILE *solution;

    char filename[255];
    char solutionfile[255];

    sprintf(solutionfile, "%s.%06d.%d", "ss.final.x", myid, myid);

    int nvalues = nx * ny;

    double sum = 0;
    double diff, diff2, L2;

    // get the localally calculated and exact solutions)
    double *values = calloc(nvalues+1, sizeof(double));

    // Opens the local solution file, reads the solution to the array, values
    // and closes the file. 

    // Opens the solution file 
    if ((solution = fopen(solutionfile, "rt")) == NULL) {
      printf("Error: can't open output file %s\n", solutionfile);
      MPI_Finalize();
      exit(1);
    }
    double mean = 0.0;
    // Reads the file to the array values 
    for (i = 0; i < nvalues + 1; i++) {
      fscanf(solution, "%lf", &values[i]);
      if (i != 0){
      mean += values[i];
      }
    }
    // Close the file 
    fflush(solution);
    fclose(solution);

    MPI_Barrier(MPI_COMM_WORLD);
mean /= nvalues; 
    // Calculates the sum of the squared differences for the exact and numerical solutions on each processor   
  for (i = 1; i < nvalues + 1; i++) {
      diff = values[i] -mean - exactsolution[i-1];
      diff2 = diff * diff;
      if (diff2 > 10){
       printf("%+6f       %+6f     %+6f\n", values[i], exactsolution[i-1],
       diff);
      }
      sum += diff2;
      }
    //printf("Myid is: %d and mysum is: %f\n", myid, sum);
    
   //  Ensure all processors are caught up 
    MPI_Barrier(MPI_COMM_WORLD);
 
    // Passes the individual sums to processor 0 and sums them. 
    double *sendbuffer, *recvbuffer;
    sendbuffer = calloc(1,  sizeof(double));
    sendbuffer[0] = sum;
    recvbuffer = calloc(1, sizeof(double));
    int root = 0;
    int count = 1;
    MPI_Reduce(sendbuffer, recvbuffer, count, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
 
  // Calculates the L2 norm  
  if (myid == 0){
      //printf("%f\n" , recvbuffer[myid]);
      // nvalues = (nx-2)*(ny-2);      
      double totalvalues = nvalues * num_procs * 1.0; 

      L2 = sqrt(recvbuffer[0] / totalvalues);      
      printf("The L2 norm is %e\n", L2);    
    }
   
    sprintf(filename, "%s.%06d", "vis/ss.sol", myid);

    if ((file = fopen(filename, "w")) == NULL) {
      printf("Error: can't open output file %s\n", filename);
      MPI_Finalize();
      exit(1);
    }

    // Save solution 
    for (i = 1; i < nvalues+1; i++) {
      fprintf(file, "%.14e\n", values[i]);
    }

    // Vis Cleanup
    fflush(file);
    fclose(file);

    free(values);
    free(exactsolution);
    free(sendbuffer);
    free(recvbuffer);

    // save global finite element mesh 
    if (myid == 0) {
      GLVis_PrintGlobalMesh("vis/ss.mesh", Nx, Ny, nx, ny, hx, hy);
    }
 } 

  // Clean up 

  HYPRE_SStructMatrixDestroy(A);
  HYPRE_SStructVectorDestroy(b);
  HYPRE_SStructVectorDestroy(x);

  // Finalize MPI
  
  MPI_Finalize();

  return (0);
}

