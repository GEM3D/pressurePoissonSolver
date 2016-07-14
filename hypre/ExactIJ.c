/* This program uses the SStruct format to set up a poisson equation, converts
 * the matrices and vectors to a sparse format and solves using AMGBOOMER
 *
 * The problem to be solved is the poisson equation where u =
 *sin(pi*x)*sin(pi*y)
 * f = -2pi^2*sin(2*pi*x)*(2*pi*y) with dirchlet boundary conditions
 *
 * Compile with: make ExactIJ
 * Execute with: mpirun -np 4 ExactIJ -nx 100 -ny 100 -NX 2 -solver 0 -vis
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
  hx = 1.0 / (Nx * nx + 1);
  hy = 1.0 / (Ny * ny + 1);

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

  // Set up the grid 
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

  // Set the dirchelet boundary condition
  {
    int bc_ilower[2];
    int bc_iupper[2];
    int nentries = 1;
    int nvaluesx = nentries * nx; /*  number of stencil entries times the length
                                      of one side of my grid box */
    int nvaluesy = nentries * ny;
    double *valuesx, *valuesy;
    int stencil_indices[1];

    valuesx = calloc(nvaluesx, sizeof(double));
    for (j = 0; j < nvaluesx; j++)
      valuesx[j] = 0.0;

    valuesy = calloc(nvaluesy, sizeof(double));
    for (j = 0; j < nvaluesy; j++)
      valuesy[j] = 0.0;

    // Recall: pi and pj describe position in the processor grid 
    if (pj == 0) {
      // Bottom row of grid points 
      bc_ilower[0] = pi * nx;
      bc_ilower[1] = pj * ny;

      bc_iupper[0] = bc_ilower[0] + nx - 1;
      bc_iupper[1] = bc_ilower[1];

      stencil_indices[0] = 3;

      HYPRE_SStructMatrixSetBoxValues(A, part, bc_ilower, bc_iupper, var,
                                      nentries, stencil_indices, valuesy);
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
      // printf("right collumn zeroed\n");
    }
    free(valuesx);
    free(valuesy);
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
    double y = (j + 1.) * 1.0 * hy + ilower[1] * hy;
    for (i = 0; i < nx; i++) {
      double x = (i + 1.) * 1.0 * hx + ilower[0] * hx;
      rhs_values[i + j * nx] =
          -8. * PI2 * sin(2. * M_PI * x) * sin(2. * M_PI * y);
      exactsolution[i + j * nx] = 1.0 * sin(2. * M_PI * x) * sin(2. * M_PI * y);
      x_values[i + j * nx] = 0.0;
    }
  }
  // printf("RHS calculation ended\n");
  HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, var, rhs_values);
  HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, var, x_values);


    printf("%d  %f\n", (ny*nx-1), exactsolution[ny*nx-1]);
  



  free(x_values);
  free(rhs_values);

  HYPRE_SStructVectorAssemble(b);
  // HYPRE_SStructVectorPrint("ss.initial.b", b, 0);
  HYPRE_SStructVectorGetObject(b, (void **)&par_b);
  HYPRE_SStructVectorAssemble(x);
  HYPRE_SStructVectorGetObject(x, (void **)&par_x);

  // Print out initial vectors 
  HYPRE_SStructMatrixPrint("SStructExact/ss.initial.A", A, 0);
  HYPRE_SStructVectorPrint("SStructExact/ss.initial.b", b, 0);
  // printf("About to solve\n");

  // Select a solver 
  if (solver_id == 0) {
    int num_iterations;
    double final_res_norm;

    // Create AMG solver 
    HYPRE_BoomerAMGCreate(&solver);

    // Set solver parameters 
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);
    HYPRE_BoomerAMGSetCoarsenType(solver, 6);
    HYPRE_BoomerAMGSetRelaxType(solver, 3);
    HYPRE_BoomerAMGSetNumSweeps(solver, 1);
    HYPRE_BoomerAMGSetMaxLevels(solver, 20);
    HYPRE_BoomerAMGSetTol(solver, 1e-12);

    // Set up and solve 
    HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
    HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);
    //printf("solution obtained\n");
  
    // Run information 
    HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    if (myid == 0) {
      printf("\n");
      printf("Iterations = %d\n", num_iterations);
      printf("Final Relative Residual Norm = %e\n", final_res_norm);
      printf("\n");
    }

    // Destroy solver 
    HYPRE_BoomerAMGDestroy(solver);
    // HYPRE_ParVectorPrint(par_x, "SStructExact/ss.final.x");
    char solutionfile[255];
    sprintf(solutionfile, "%s.%06d", "ss.final.x", myid);
    HYPRE_ParVectorPrint(par_x, solutionfile);
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

    // Reads the file to the array values 
    for (i = 0; i < nvalues + 1; i++) {
      fscanf(solution, "%lf", &values[i]);
    }

    // Close the file 
    fflush(solution);
    fclose(solution);

    MPI_Barrier(MPI_COMM_WORLD);

    // Calcluates the sum of the squared differences for the exact and numerical solutions on each processor   
  for (i = 1; i < nvalues + 1; i++) {
      diff = values[i] - exactsolution[i-1];
      diff2 = diff * diff;
      if (diff2 > 10){
       printf("%+6f       %+6f     %+6d\n", values[i], exactsolution[i-1],
       i-1);
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
    //  nvalues = (nx-2)*(ny-2);      
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

  // Finalize MPI*/
  
  MPI_Finalize();

  return (0);
}

