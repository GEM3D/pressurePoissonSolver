/* This program uses the SStruct format to set up a poisson equation on a semi-structured mesh.
 *
 * The poisson equation solved is:
 * u = sin(pi*x)*sin(pi*y)
 * where
 * f = -2pi^2*sin(2*pi*x)*sin(2*pi*y)
 * with neumann boundary conditions
 *
 * compile with: make two_parts_ss
 * Execute with: mpirun -np 1 two_parts_ss
 *
 * Created by: Joseph McNeal
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_parcsr_ls.h"
#include "vis.c"
#include "add_to_vis.c"
#include "HYPRE_parcsr_mv.h"
#include "mpi.h"
#define  PI2 M_PI*M_PI

int main(int argc, char *argv[]){
  //printf("start\n");

  //Initialize Variables
  int i, j, pi, pj;
  int myid, num_procs, num_parts;
  int Nx, Ny; //Nx and Ny are the number of Processors in each direction
  int nx, ny; //nx and ny are number of mesh cells per proc
  int solver_id;
  int vis;
  int ilower[2], iupper[2];
  double hx, hy;
  //hx and hy are the default grid spaces
  //phx and phy are the grid spacings for each part
 
  //printf("c variables initialized\n");

  HYPRE_SStructGraph graph;
  HYPRE_SStructGrid grid;
  HYPRE_SStructStencil stencil;
  HYPRE_SStructMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_SStructVector b;
  HYPRE_ParVector par_b;
  HYPRE_SStructVector x;
  HYPRE_ParVector par_x;

  //printf("hypre variables initialized\n");
  
  //Set problem parameters
  Nx = 1;
  Ny = 1;  //Change this when more than one processor needs to be used
  num_parts = 2;
  nx = 10;
  ny = 10;
  solver_id = 0;
  vis = 1;

  int plower[2*num_parts], pupper[2*num_parts];
  //int plowerloc[2],pupperloc[2];
  int part[num_parts];
  double phx[num_parts], phy[num_parts];

  for(i = 0; i < num_parts; i++){
      part[i] = i;
  }

  //printf("variables set\n");

   //Initialize MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  //printf("mpi initialized\n");


  //Set up processors in grid 
  //Determine the global grid spacing
  hx = 1.0 / 10;
  hy = 1.0 / 10;
  

  //Determine the grid spacing of each part
  //Here Part 0 has the global grid spacing, and part 1 has half of that
  //Set the parts into a grid
  int num_partsx = num_parts;
  int num_partsy = 1;
  phx[0] = hx;
  phy[0] = hy;
  phx[1] = hx;
  phy[1] = hy;

  if (Ny <= Nx) {
       pj = myid / Nx;
       pi = myid - pj * Ny;
       } else {
       pi = myid / Ny;
       pj = myid - pi * Nx;
       }

   //printf("processors assigned\n");

  // Determine each processor's piece of the grid 
  //ilower[0] = pi * nx;
  //ilower[1] = pj * ny;

  //iupper[0] = ilower[0] + nx - 1;
  //iupper[1] = ilower[1] + ny - 1;
       ilower[0] = 1;
       ilower[1] = 1;
       iupper[0] = 10;
       iupper[1] = 10;
  
  //printf("grid divided to processors\n");

   HYPRE_SStructGridCreate(MPI_COMM_WORLD, 2, num_parts, &grid);

  //Determine each part's piece of the processor's grid
  

  for(i = 0; i < num_parts; i++){
       plower[2 * i] = ilower[0] + i * nx/num_partsx;
       plower[1 + 2 * i] = ilower[1];
       //plowerloc[0] = plower[2*i];
       //plowerloc[1] = plower[1+2*i];
       pupper[2 * i] = plower[0 + 2 * i] + nx/num_partsx - 1;
       pupper[1 + 2 * i] = plower[1 + 2 * i] + ny/num_partsy - 1;
       //pupperloc[0] = pupper[2*i];
       //pupperloc[1] = pupper[1+2*i];

       //HYPRE_SStructGridSetExtents(grid, part[i], plowerloc, pupperloc);
   }

   plower[0] = 1;
   plower[1] = 1; 
   pupper[0] = 10;
   pupper[1] = 10;
   plower[2] = 11;
   plower[3] = 1;
   pupper[2] = 20;
   pupper[3] = 10;

   HYPRE_SStructGridSetExtents(grid, 0, ilower, iupper);
   ilower[0] = 11;
   ilower[1] = 1;
   iupper[0] = 20;
   iupper[1] = 10;

   HYPRE_SStructGridSetExtents(grid, 1, ilower, iupper);
   //printf("The coordinates for the grid are [%d, %d] [%d, %d]\n", ilower[0], ilower[1], iupper[0], iupper[1]);
   printf("The coordinates for part 0 are [%d, %d] [%d, %d]\n", plower[0], plower[1], pupper[0], pupper[1]);
   printf("The coordinates for part 1 are [%d, %d] [%d, %d]\n", plower[2], plower[3], pupper[2], pupper[3]);

   //Set variable Conditions
   int nvars = 1;
   HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
   //printf("Variable type set\n");
 
   for(i = 0; i < num_parts; i++){
       HYPRE_SStructGridSetVariables(grid, part[i] , nvars, vartypes);
   }
 
   //printf("part variables set\n");
 

   //Set the relationship between the parts
 
   i = 0;
   int nbor_part = 1;
   int b_ilower[2] = {11, 1};
   int b_iupper[2] = {11, 10};
   int indexmap[2] = {0,1};
   int index_dir[2] = {1,1};

   HYPRE_SStructGridSetNeighborPart(grid, i, b_ilower, b_iupper, nbor_part, b_ilower, b_iupper, indexmap, index_dir);

   i = 1;
   nbor_part = 0;

   b_ilower[0] = 10;
   b_ilower[1] = 1;
   b_iupper[0] = 10;
   b_iupper[1] =  10;
   int indexmap2[2] = {0,1};
   int index_dir2[2] = {1,1};

   HYPRE_SStructGridSetNeighborPart(grid, i, b_ilower, b_iupper, nbor_part, b_ilower, b_iupper, indexmap2, index_dir2);

   HYPRE_SStructGridAssemble(grid);

   //printf("Grid assembled\n");
 
   //Set up the Stencil
   HYPRE_SStructStencilCreate(2, 5, &stencil);
   int entry;
   int offsets[5][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1} };
   int var = 0;
 
   for (entry = 0; entry < 5; entry++){
        HYPRE_SStructStencilSetEntry(stencil, entry, offsets[entry], var);
   }

   //printf("Stencil Setup\n");

   //Create the Graph
   HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);
 
   for (i = 0; i < num_parts; i++){
      HYPRE_SStructGraphSetStencil(graph, part[i], var, stencil);
   }
 
   HYPRE_SStructGraphAssemble(graph);
   //printf("Graph Assembled\n");
   //Setup the matrix

   HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);
   HYPRE_SStructMatrixSetObjectType(A, HYPRE_PARCSR);
   HYPRE_SStructMatrixInitialize(A);
 
   int nentries = 5;
   double *values;

 
   int stencil_indices[5] = {0,1,2,3,4};
   int nvalues[num_parts];
   double phx2[num_parts], phy2[num_parts];
   double phx2inv[num_parts], phy2inv[num_parts];
   int lower[2], upper[2];
   int nvaluesx, nvaluesy, stencil_index[1], center_index[1];
   center_index[0] = 0;
   double *valuesx, *valuesy, *center_valuesx, *center_valuesy;

   //printf("more variables setup\n");

   for(i = 0; i < num_parts; i++){
       nentries = 5;
       //printf("phy[%d] = %f\n", i, phy[i]);
       nvalues[i] = (1/phx[i] * 1/phy[i])*nentries;
       //printf("entries:%d\n", nentries);
       //printf("nvalues set\n");
       phx2[i] = phx[i] * phx[i];
       //printf("phx2 set\n");
       phy2[i] = phy[i] * phy[i];
       //printf("phy2 set\n");
       phx2inv[i] = 1./phx2[i];
       //printf("phx2inv set\n");
       phy2inv[i] = 1./phy2[i];
       //printf("phy2inv set\n");
       values = calloc(nvalues[i], sizeof(double));

       MPI_Barrier(MPI_COMM_WORLD);


       //printf("%d\n",nvalues[i]);
       //printf("%f\n",phx2inv[i]);
       //printf("%f\n",phy2inv[i]);
       //printf("%f\n",values[nvalues[i]]);
 
       //Account for variation in hx and hy
       //change to weights
       for(j = 0; j < nvalues[i]-nentries; j += nentries) {
             values[j] = -2.0 * phx2inv[i] - 2.0 * phy2inv[i];
             values[j+1] = phx2inv[i];
             values[j+2] = phx2inv[i];
             values[j+3] = phy2inv[i];
             values[j+4] = phy2inv[i];
             //printf("p: %d, [%+.6f,%+.6f,%+.6f,%+.6f,%+.6f]\n", i, values[j], values[j+1], values[j+2], values[j+3], values[j+4]);
       }
       MPI_Barrier(MPI_COMM_WORLD);
       printf("%d\n", j+4);
       if (i == 0){       
             lower[0] = 1; 
             lower[1] = 1;
             //printf("lowers set\n");
             upper[0] = 10;
             upper[1] = 10;
       }

       if (i == 1){
             lower[0] = 11;
             lower[1] = 1;
             upper[0] = 20;
             upper[1] = 10;
       }
       //printf("uppers set\n");
       HYPRE_SStructMatrixSetBoxValues(A, part[i], lower, upper, var, nentries, stencil_indices, values);
       //printf("box values set\n");
       //Set the Neumann Boundary Conditions
       //Check for different numbers of parts
       nentries = 1;
       nvaluesx = nentries / phx[i];
       nvaluesy = nentries / phy[i];
       printf("nx: %d, ny: %d\n", nvaluesx, nvaluesy); 
       valuesx = calloc(nvaluesx, sizeof(double));
       valuesy = calloc(nvaluesy, sizeof(double));
       center_valuesx = calloc(nvaluesx, sizeof(double));
       center_valuesy = calloc(nvaluesy, sizeof(double));
       for(j = 0; j < nvaluesx; j++){
              valuesx[j] = 0.0;
              center_valuesx[j] = phx2inv[i];
       }
 
       for(j = 0; j < nvaluesy; j++){
              valuesy[j] = 0.0;
              center_valuesy[j] = phy2inv[i];
       }

       //printf("center values set\n");

       if (pj == 0){
              //Bottom Row of grid points
              b_ilower[0] = plower[2*i];
              b_ilower[1] = plower[1+2*i];
              b_iupper[0] = pupper[1+2*i];
              b_iupper[1] = b_ilower[1];
              stencil_index[0] = 3;
              if (myid == 0 && i == 0 ){
                     center_valuesy[1] = -phy2inv[0];
              }
              HYPRE_SStructMatrixSetBoxValues(A, part[i], b_ilower, b_iupper, var, nentries, stencil_index, valuesy);
              HYPRE_SStructMatrixAddToBoxValues(A, part[i], b_ilower, b_ilower, var, nentries, center_index, center_valuesy);

              center_valuesy[1] = phy2inv[i];
       }


       if (pj == Ny - 1) {
              //Upper row of grid points
              b_ilower[0] = plower[2*i];
              b_ilower[1] = pupper[1+2*i];

              b_iupper[0] = pupper[2*i];
              b_iupper[1] = b_ilower[1];

              stencil_index[0] = 4;

              HYPRE_SStructMatrixSetBoxValues(A, part[i], b_ilower, b_iupper, var,
                                      nentries, stencil_index, valuesy);
              HYPRE_SStructMatrixAddToBoxValues(A, part[i], b_ilower, b_iupper, var, nentries, center_index, center_valuesy);
       }

       if (i == 0) {
               //printf("foo\n");
               //Left column of gridpoints
               b_ilower[0] = plower[2*i];
               b_ilower[1] = plower[1+2*i];

               b_iupper[0] = b_ilower[0];
               b_iupper[1] = pupper[1+2*i];

               stencil_index[0] = 1;
       
                HYPRE_SStructMatrixSetBoxValues(A, part[i], b_ilower, b_iupper, var, nentries, stencil_index, valuesx);
                HYPRE_SStructMatrixAddToBoxValues(A, part[i], b_ilower, b_iupper, var, nentries, center_index, center_valuesx);
       }

       if (i == 1) {
               //Right column of gridpoints
                b_ilower[0] = pupper[2*i];
                b_ilower[1] = plower[1+2*i];

                b_iupper[0] = b_ilower[0];
                b_iupper[1] = pupper[1+2*i];

                stencil_index[0] = 2;

                HYPRE_SStructMatrixSetBoxValues(A, part[i], b_ilower, b_iupper, var, nentries, stencil_index, valuesx);
                HYPRE_SStructMatrixAddToBoxValues(A, part[i], b_ilower, b_iupper, var, nentries, center_index, center_valuesx);
       }
 
       free(valuesx);
       free(center_valuesx);
       free(valuesy);
       free(center_valuesy);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  //printf("BC set\n");
  HYPRE_SStructMatrixAssemble(A);
  HYPRE_SStructMatrixGetObject(A, (void **)&parcsr_A);

  //Create the RHS and Solution Vectors
  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &b);
  HYPRE_SStructVectorSetObjectType(b, HYPRE_PARCSR);
  HYPRE_SStructVectorInitialize(b);

  HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &x);
  HYPRE_SStructVectorSetObjectType(x, HYPRE_PARCSR);
  HYPRE_SStructVectorInitialize(x);

  
  //Set the RHS, Exact Solution, and initial guess values
  double *rhs_values, *x_values, *exactsolution;
  int numvalues = 0;
  int *part_x, *part_y;
  part_x = calloc(num_parts, sizeof(int));
  part_y = calloc(num_parts, sizeof(int));
  for(i = 0; i < num_parts; i++) {
      //numvalues += (pupper[2*i] - plower[2*i]+1)/phx[i];
      part_x[i] = (pupper[2*i] - plower[2*i]+1)*(hx/phx[i]); //look at rounding
      //numvalues += (pupper[1+2*i] - plower[1+2*i]+1)/phy[i];
      part_y[i] += (pupper[1+2*i] - plower[1+2*i]+1)*(hy/phy[i]);
        numvalues += part_x[i]*part_y[i];
      //printf("(%d - %d)/%f\n", pupper[1+2*i], plower[1+2*i], phy[i]);
      //printf("%d\n", part_x[i]);
      //printf("%d\n", part_y[i]);
  }
 
  MPI_Barrier(MPI_COMM_WORLD);
  //printf("%d\n", numvalues);
  //printf("%d + %d\n", part_x[0], part_y[0]);
  //printf("%d + %d\n", part_x[1], part_y[1]);
  rhs_values = calloc(nx*ny, sizeof(double));//Change later
  x_values = calloc(nx*ny, sizeof(double));
  exactsolution = calloc(numvalues, sizeof(double));

  //printf("Vector memory allocated\n");
  //printf("(%d,%d),(%d,%d)\n", plower[0], plower[1], pupper[0], pupper[1]);
  //printf("(%d,%d),(%d,%d)\n", plower[2], plower[3], pupper[2], pupper[3]);

  int p;
  int counter = 0;
  for(p = 0; p < num_parts; p++){

        //Get the boundaries for the part
        if (p == 0){
           lower[0] = 1;
           lower[1] = 1;
           upper[0] = nx;
           upper[1] = ny;
        }

        if (p == 1){
           lower[0] += nx;
           upper[0] += nx;
        }
        printf("Lowers: [%d, %d]\n Uppers: [%d,%d]\n",lower[0], lower[1], upper[0], upper[1]);
        printf("pupper[%d]: %d\n",p, pupper[2*p]);
        printf("hx: %f,  hy: %f\n", phx[p], phy[p]);
                        
  for(j = 0; j < ny; j++){
        double y = (j + .5) * hy;
        for(i = 0; i < nx; i++){
               double x = p + (i +.5) * hx;
               rhs_values[i + j * nx] = -8.* PI2 * cos(2. * M_PI * x) * cos(2. * M_PI * y);
               exactsolution[i + j * nx] = 1.0 * cos(2. * M_PI * x) * cos(2.* M_PI * y);
               x_values[i + j * nx] = 1.0;
               counter++;
               //printf("rhs: %f\n", rhs_values[i + j * nx]);
               printf("p: %d, x: %f, y: %f, rhs: %+.6f  exact: %+.6f\n", p, x, y, rhs_values[i+j*nx], exactsolution[i+j*nx]);
         }
  
        //printf("total number of points:%d\n", counter);
     }


        //Set the rhs and initial guess values into the hypre vectors
        HYPRE_SStructVectorSetBoxValues(b, p, lower, upper, var, rhs_values);
        //printf("p: %d  l: [%d,%d]  u: [%d,%d]\n", part[p], lower[0], lower[1], upper[0], upper[1]);
        //HYPRE_SStructVectorPrint("multipart_data/ss.wut.b", b, 0);
        HYPRE_SStructVectorSetBoxValues(x, p, lower, upper, var, x_values);
        //printf("vectors actually set\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
  free(x_values);
  free(rhs_values);

  //printf("temp vectors freed\n");

  HYPRE_SStructVectorAssemble(b);
  HYPRE_SStructVectorGetObject(b, (void **)&par_b);
  HYPRE_SStructVectorAssemble(x);
  HYPRE_SStructVectorGetObject(x, (void **)&par_x);
  //printf("vectors assembled\n");

  //Print out the RHS vector and A matrix
  HYPRE_SStructMatrixPrint("multipart_data/ss.initial.A", A, 0);
  HYPRE_SStructVectorPrint("multipart_data/ss.initial.b", b, 0);
  
  //printf("vectors printed out\n");

  int num_iterations;
  double final_res_norm, t1, t2;

  //Select a solver
  if (solver_id == 0){
         //Create AMG solver
         HYPRE_Solver solver;
         HYPRE_BoomerAMGCreate(&solver);
         
         //Set solver parameters
         HYPRE_BoomerAMGSetPrintLevel(solver, 3);
         HYPRE_BoomerAMGSetStrongThreshold(solver, .15);
         HYPRE_BoomerAMGSetCoarsenType(solver, 6);
         HYPRE_BoomerAMGSetRelaxType(solver, 6);
         HYPRE_BoomerAMGSetNumSweeps(solver, 1);
         HYPRE_BoomerAMGSetMaxLevels(solver, 20);
         HYPRE_BoomerAMGSetTol(solver, 1e-12);
         HYPRE_BoomerAMGSetInterpType(solver, 14);
         HYPRE_BoomerAMGSetAggInterpType(solver, 2);
         HYPRE_BoomerAMGSetMaxIter(solver, 100);
         //printf("solver parameters set\n");
         HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
         //printf("solver setup\n");
         t1 = MPI_Wtime();
         HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);
         t2 = MPI_Wtime();
         //printf("solve complete\n");
   
         //Get Run Information
         HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
         HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

         HYPRE_BoomerAMGDestroy(solver);
   }



   if (solver_id == 1) {

       //Set the PCG solvers parameters  
       HYPRE_Solver cgsolver;
       HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &cgsolver);
       HYPRE_ParCSRPCGSetTol(cgsolver, 1.0e-12);
       HYPRE_ParCSRPCGSetMaxIter(cgsolver, 500);
       HYPRE_ParCSRPCGSetTwoNorm(cgsolver, 1);
       HYPRE_ParCSRPCGSetLogging(cgsolver, 3);

       //Set the AMG preconditioner parameters  
       HYPRE_Solver precond;
       HYPRE_BoomerAMGCreate(&precond);
       HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
       HYPRE_BoomerAMGSetCoarsenType(precond, 6);
       HYPRE_BoomerAMGSetTol(precond, 0.0);
       //HYPRE_BoomerAMGSetPrintLevel(precond, 1);
       HYPRE_BoomerAMGSetMaxIter(precond, 1);

       // Setup and Solve the system
       HYPRE_ParCSRPCGSetPrecond(cgsolver, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, precond);
       HYPRE_ParCSRPCGSetup(cgsolver, parcsr_A, par_b, par_x);
       t1 = MPI_Wtime();
       HYPRE_ParCSRPCGSolve(cgsolver, parcsr_A, par_b, par_x);
       t2 = MPI_Wtime();

       //Get run info
       HYPRE_ParCSRPCGGetNumIterations(cgsolver, &num_iterations);
       HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(cgsolver, &final_res_norm);

       //Clean up
       HYPRE_BoomerAMGDestroy(precond);
       HYPRE_ParCSRPCGDestroy(cgsolver);
   }


   if (myid == 0) {
           printf("\n");
           printf("Iterations = %d\n", num_iterations);
           printf("Final Relative Residual Norm = %e\n", final_res_norm);
           printf("\n");
   }


   MPI_Barrier(MPI_COMM_WORLD);

   double *time, *timemax;
   time = calloc(1,  sizeof(double));
   time[0] = t2-t1;
   timemax = calloc(1, sizeof(double));
   int root = 0;
   int count = 1;
   MPI_Reduce(time, timemax, count, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);

   if (myid == 0){
           printf("Elapsed Wall Time %f\n", timemax[0]);
   }


   char solutionfile[255];
   sprintf(solutionfile, "%s.%06d", "mp.solution.x", myid);
   HYPRE_ParVectorPrint(par_x, solutionfile);


   //Save the solution for visualization and L2 Norm Calculation
   if (vis) {
        FILE *file;
        FILE *solution;

        char filename[255];
        char solutionfile[255];

        sprintf(solutionfile, "%s.%06d.%d", "mp.solution.x", myid, myid);

        int numvalues = 0;
        int xval = 0;
        int yval = 0;
        int partval = 0;
        for(i = 0; i < num_parts; i++) {
             //Change this if size isn't 1
             xval = 1/phx[i];
             yval = 1/phy[i];
             partval = xval * yval;
             //printf("partnum: %d\n", partval);
             numvalues += partval;
        }


        double sum = 0;
        double diff, diff2, L2;
 
        //Get the localally calculated and exact solutions
        double *values = calloc(numvalues+1, sizeof(double));

        //Opens the local solution file, reads the solution to the array values and then closes the file

        //Opens the solution file
        if ((solution = fopen(solutionfile, "rt")) == NULL) {
             printf("Error: can't open output file %s\n", solutionfile);
             MPI_Finalize();
             exit(1);
        }

        double mean = 0.0;
        //Reads the file to the array values
        for (i = 0; i < numvalues + 1; i++) {
              fscanf(solution, "%lf", &values[i]);
              if (i != 0){
                     mean += values[i];
              }
        }
        printf("mean: %f\n", mean);
        //Close the file
        fflush(solution);
        fclose(solution);

        MPI_Barrier(MPI_COMM_WORLD);
        double *sendmean, *recvmean;
        recvmean = calloc(1, sizeof(double));
        sendmean = calloc(1, sizeof(double));
        sendmean[0] = mean;
        MPI_Allreduce(sendmean, recvmean, count ,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        recvmean[0] /= (numvalues*num_procs);

        //Calculates the sum of the squared differences for the exact and numerical solution on each processor
        for (i = 1; i < numvalues + 1; i++) {
             diff = values[i] -recvmean[0] - exactsolution[i-1];
             diff2 = diff * diff;
             sum += diff2;
        }
       
        MPI_Barrier(MPI_COMM_WORLD);
    
       //Passes the individual sums to processor 0 and sums them
       double *sendbuffer, *recvbuffer;
       sendbuffer = calloc(1,  sizeof(double));
       sendbuffer[0] = sum;
       recvbuffer = calloc(1, sizeof(double));
       int root = 0;
       int count = 1;
       MPI_Reduce(sendbuffer, recvbuffer, count, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
       MPI_Barrier(MPI_COMM_WORLD);

       //Calculates the L2 norm
       if (myid == 0){
           double totalvalues = numvalues * num_procs * 1.0;
           printf("L2 values: %f\n", totalvalues);

           L2 = sqrt(recvbuffer[0] / totalvalues);
           printf("The L2 norm is %e\n", L2);
       }

       sprintf(filename, "%s.%06d", "multipart_data/mp.sol", myid);

       if ((file = fopen(filename, "w")) == NULL) {
           printf("Error: can't open output file %s\n", filename);
           MPI_Finalize();
           exit(1);
       }

       //Save the solution
       for (i = 1; i < numvalues+1; i++) {
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
       if (myid == 0) {
       //GLVis_PrintGlobalMesh("multipart_data/ss.mesh", Nx, Ny, nx, ny, hx, hy);
       }
   }



   //Cleanup
   HYPRE_SStructMatrixDestroy(A);
   HYPRE_SStructVectorDestroy(b);
   HYPRE_SStructVectorDestroy(x);
  
   MPI_Finalize();

   return(0);
}









  
