/*

   Interface:      Structured interface (Struct)

   Compile with:   make try2varxy

   Execute with:   mpirun -np 1 try2_varxy -nx 101 -ny 101 -solver 0 -v 1 1 -vis -Nx 1



*/

#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
#include "add_to_vis.c"
#include "vis.c"
#define PI2 M_PI*M_PI

int main (int argc, char *argv[])
{
   int i, j, k;

   int myid, num_procs;

   int nx, ny, Nx, Ny, pi, pj;
   double hx, hy;
   int ilower[2], iupper[2];

   int solver_id;
   int n_pre, n_post;

   HYPRE_StructGrid     grid;
   HYPRE_StructStencil  stencil;
   HYPRE_StructMatrix   A;
   HYPRE_StructVector   b;
   HYPRE_StructVector   analytical;
   HYPRE_StructVector   x;
   HYPRE_StructSolver   solver;
   HYPRE_StructSolver   precond;
   
   int num_iterations;
   double final_res_norm;

   int vis;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   /* Set defaults */
   nx = 33;
   ny = 33;
   Nx = 1;
   solver_id = 0;
   n_pre  = 1;
   n_post = 1;
   vis = 0;

   /* Parse command line */
   {
      int arg_index = 0;
      int print_usage = 0;

      while (arg_index < argc)
      {
         if ( strcmp(argv[arg_index], "-nx") == 0 )
         {
            arg_index++;
            nx = atoi(argv[arg_index++]);
         }

         else if ( strcmp(argv[arg_index], "-ny") == 0)
         {
            arg_index++;
            ny = atoi(argv[arg_index++]);
         }


         else if ( strcmp(argv[arg_index], "-Nx") == 0)
         {
            arg_index++;
            Nx = atoi(argv[arg_index++]);
         }

         else if ( strcmp(argv[arg_index], "-solver") == 0 )
         {
            arg_index++;
            solver_id = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-v") == 0 )
         {
            arg_index++;
            n_pre = atoi(argv[arg_index++]);
            n_post = atoi(argv[arg_index++]);
         }
         else if ( strcmp(argv[arg_index], "-vis") == 0 )
         {
            arg_index++;
            vis = 1;
         }
         else if ( strcmp(argv[arg_index], "-help") == 0 )
         {
            print_usage = 1;
            break;
         }
         else
         {
            arg_index++;
         }
      }

      if ((print_usage) && (myid == 0))
      {
         printf("\n");
         printf("Usage: %s [<options>]\n", argv[0]);
         printf("\n");
         printf("  -n <n>              : problem size per processor (default: 33)\n");
         printf("  -solver <ID>        : solver ID\n");
         printf("                        0  - PCG with SMG precond (default)\n");
         printf("                        1  - SMG\n");
         printf("  -v <n_pre> <n_post> : number of pre and post relaxations (default: 1 1)\n");
         printf("  -vis                : save the solution for GLVis visualization\n");
         printf("\n");
      }

      if (print_usage)
      {
         MPI_Finalize();
         return (0);
      }
   }

   /* Figure out the processor grid (N x N).  The local problem
      size for the interior nodes is indicated by n (n x n).
      pi and pj indicate position in the processor grid. */
   Ny  = num_procs/Nx;
   hx  = 1.0 / (Nx*nx+1); /* note that when calculating h we must
                          remember to count the boundary nodes */
   hy  = 1.0/ (Ny*ny+1);
   if (Ny <= Nx)
   {
        pj = myid / Nx;
        pi = myid - pj*Ny;
   }
   else
   {
        pi = myid / Ny;
        pj = myid - pi*Nx;
   }

   printf("MyId: %d, np: %d, Nx: %d, Ny: %d, hx: %f, hy: %f, pj: %d, pi: %d\n",myid, num_procs, Nx, Ny, hx, hy, pj, pi);
   /* Figure out the extents of each processor's piece of the grid. */
   ilower[0] = pi*nx;
   ilower[1] = pj*ny;

   iupper[0] = ilower[0] + nx-1;
   iupper[1] = ilower[1] + ny-1;

   /* 1. Set up a grid */
   {
      /* Create an empty 2D grid object */
      HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);

      /* Add a new box to the grid */
      HYPRE_StructGridSetExtents(grid, ilower, iupper);

      /* This is a collective call finalizing the grid assembly.
         The grid is now ``ready to be used'' */
      HYPRE_StructGridAssemble(grid);
   }

   /* 2. Define the discretization stencil */
   {
      /* Create an empty 2D, 5-pt stencil object */
      HYPRE_StructStencilCreate(2, 5, &stencil);

      /* Define the geometry of the stencil */
      {
         int entry;
         int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};

         for (entry = 0; entry < 5; entry++)
            HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
      }
   }

   /* 3. Set up a Struct Matrix */
   {
      int nentries = 5;
      int nvalues = nentries*nx*ny;
      //printf("%d\n", nvalues);
      double *values;
      int stencil_indices[5];

      /* Create an empty matrix object */
      HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);

      /* Indicate that the matrix coefficients are ready to be set */
      HYPRE_StructMatrixInitialize(A);

      values = calloc(nvalues, sizeof(double));

      for (j = 0; j < nentries; j++)
         stencil_indices[j] = j;

      /* Set the standard stencil at each grid point,
         we will fix the boundaries later */

     double  hx2 = hx * hx;
     double  hy2 = hy * hy;
     double  hx2inv = 1/hx2;
     double  hy2inv = 1/hy2;

      for (i = 0; i < nvalues; i += nentries)
      {
         values[i] = -2.0*hx2inv-2.0*hy2inv;
         values[i+1] = hx2inv;
         values[i+2] = hx2inv;
         values[i+3] = hy2inv;
         values[i+4] = hy2inv;
      }

      HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                     stencil_indices, values);

      free(values);
   }

   /* 4. Incorporate the zero boundary conditions: go along each edge of
         the domain and set the stencil entry that reaches to the boundary to
         zero.*/
   {
      int bc_ilower[2];
      int bc_iupper[2];
      int nentries = 1;
      int nvaluesx  = nentries*nx; /*  number of stencil entries times the length     
                               of one side of my grid box */
      int nvaluesy  = nentries*ny;
      double *valuesx, *valuesy;
      int stencil_indices[1];

      valuesx = calloc(nvaluesx, sizeof(double));
      for (j = 0; j < nvaluesx; j++)
         valuesx[j] = 0.0;
      
      valuesy = calloc(nvaluesy, sizeof(double));
      for (j =0; j < nvaluesy; j++)
         valuesy[j] = 0.0;

      /* Recall: pi and pj describe position in the processor grid */
      if (pj == 0)
      {
         /* Bottom row of grid points */
         bc_ilower[0] = pi*nx;
         bc_ilower[1] = pj*ny;

         bc_iupper[0] = bc_ilower[0] + nx-1;
         bc_iupper[1] = bc_ilower[1];

         stencil_indices[0] = 3;

         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, valuesy);
         //printf("Bottom row zeroed\n");
      }

      if (pj == Ny-1)
      {
         /* upper row of grid points */
         bc_ilower[0] = pi*nx;
         bc_ilower[1] = pj*ny + ny-1;

         bc_iupper[0] = bc_ilower[0] + ny-1;
         bc_iupper[1] = bc_ilower[1];

         stencil_indices[0] = 4;

         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, valuesy);
         //printf("Upper row zeroed\n");
      }

      if (pi == 0)
      {
         /* Left row of grid points */
         bc_ilower[0] = pi*nx;
         bc_ilower[1] = pj*ny;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + ny-1;

         stencil_indices[0] = 1;

         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, valuesx);
         //printf("left collumn zeroed\n");
      }

      if (pi == Nx-1)
      {
         /* Right row of grid points */
         bc_ilower[0] = pi*nx + nx-1;
         bc_ilower[1] = pj*ny;

         bc_iupper[0] = bc_ilower[0];
         bc_iupper[1] = bc_ilower[1] + ny-1;

         stencil_indices[0] = 2;

         HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                        stencil_indices, valuesx);
         //printf("right collumn zeroed\n");
      }

      free(valuesx);
      free(valuesy);
   }

   /* This is a collective call finalizing the matrix assembly.
      The matrix is now ``ready to be used'' */
   HYPRE_StructMatrixAssemble(A);

   /* 5. Set up Struct Vectors for b and x */
   {
      int    nvalues = nx*ny;
      double *values, *exact_solution,*xvalues;

      values = calloc(nvalues, sizeof(double));
      exact_solution = calloc(nvalues, sizeof(double));
      xvalues = calloc(nvalues, sizeof(double));

      /* Create an empty vector object */
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
      HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &analytical);

      /* Indicate that the vector coefficients are ready to be set */
      HYPRE_StructVectorInitialize(b);
      HYPRE_StructVectorInitialize(x);
      HYPRE_StructVectorInitialize(analytical);
    
     //printf("My Id is %d. %d,%d Lowers are %d,  %d and uppers are %d, %d\n",myid, pi, pj, ilower[0], ilower[1] , iupper[0], iupper[1]);
/*
     Set the values 
      for (j = ilower[1]; j <= iupper[1]; j++){
          double y = (j+1.) * 1.0 * hy;
          for (i = ilower[0]; i <= iupper[0]; i ++){
              double x = (i+1.) * 1.0 * hx;
              exact_solution[i + j * nx] = 1.0 * sin(2 * M_PI * x) * sin(2 * M_PI * y);
              values[i + j * nx] = -8. * PI2 * sin(2. * M_PI * x) * sin(2. * M_PI * y);
              //printf("my id is %d (x,y) %f, %f\n",myid, x, y);
              xvalues[i + j * nx] = 0.0;
              printf("At %f, %f the value is %f and exact solution is %f\n", x , y , values[i+j*nx], exact_solution[i+j*nx]);
          }
      }
*/

   for (j = 0; j < ny; j++){
          double y = (j+1.) * 1.0 * hy + ilower[1] * hy;
          for (i = 0; i < nx; i ++){
              double x = (i+1.) * 1.0 * hx + ilower[0] * hx;
              exact_solution[i + j * nx] = 1.0 * sin(2 * M_PI * x) * sin(2 * M_PI * y);
              values[i + j * nx] = -8. * PI2 * sin(2. * M_PI * x) * sin(2. * M_PI * y);
              //printf("my id is %d (x,y) %f, %f\n",myid, x, y);
              xvalues[i + j * nx] = 0.0;
              // printf("At %f, %f the value is %f and exact solution is %f\n", x , y , values[i+j*nx], exact_solution[i+j*nx]);
                                                    }
                                                          }
              
      HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);
      HYPRE_StructVectorSetBoxValues(analytical, ilower, iupper, exact_solution);
      HYPRE_StructVectorSetBoxValues(x, ilower, iupper, xvalues);

      free(values);
      free(xvalues);


      /* This is a collective call finalizing the vector assembly.
         The vector is now ``ready to be used'' */
      HYPRE_StructVectorAssemble(b);
      HYPRE_StructVectorAssemble(analytical);
      HYPRE_StructVectorAssemble(x);
      HYPRE_StructVectorPrint("try2varxy/try2varxy.vect.b", b, 0);
      HYPRE_StructVectorPrint("try2varxy/try2varxy.vect.exact", analytical, 0);
      HYPRE_StructVectorPrint("try2varxy/try2varxy.vect.x", x, 0);
      HYPRE_StructMatrixPrint("try2varxy/try2varxy.matrix.A", A, 0);
   }
   printf("vectors assembled\n");   

   /* 6. Set up and use a struct solver
      (Solver options can be found in the Reference Manual.) */
   if (solver_id == 0)
   {
      HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
      HYPRE_StructPCGSetMaxIter(solver, 50 );
      HYPRE_StructPCGSetTol(solver, 1.0e-12 );
      HYPRE_StructPCGSetTwoNorm(solver, 1 );
      HYPRE_StructPCGSetRelChange(solver, 0 );
      HYPRE_StructPCGSetPrintLevel(solver, 2 ); /* print each CG iteration */
      HYPRE_StructPCGSetLogging(solver, 1);
      
      //printf("solver setup\n");

      /* Use symmetric SMG as preconditioner */
      HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
      HYPRE_StructSMGSetMemoryUse(precond, 0);
      HYPRE_StructSMGSetMaxIter(precond, 1);
      HYPRE_StructSMGSetTol(precond, 0.0);
      HYPRE_StructSMGSetZeroGuess(precond);
      HYPRE_StructSMGSetNumPreRelax(precond, 1);
      HYPRE_StructSMGSetNumPostRelax(precond, 1);
      
      //printf("preconditioner setup\n");

      /* Set the preconditioner and solve */
      HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve,
                                  HYPRE_StructSMGSetup, precond);
    
      //printf("precond set\n");  

      HYPRE_StructPCGSetup(solver, A, b, x);

      //printf("pcg setup\n");

      HYPRE_StructPCGSolve(solver, A, b, x);
      
      //printf("solved\n");

      /* Get some info on the run */
      HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
      HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      /* Clean up */
      HYPRE_StructPCGDestroy(solver);
   }

   if (solver_id == 1)
   {
      HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
      HYPRE_StructSMGSetMemoryUse(solver, 0);
      HYPRE_StructSMGSetMaxIter(solver, 50);
      HYPRE_StructSMGSetTol(solver, 1.0e-12);
      HYPRE_StructSMGSetRelChange(solver, 0);
      HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
      HYPRE_StructSMGSetNumPostRelax(solver, n_post);
      /* Logging must be on to get iterations and residual norm info below */
      HYPRE_StructSMGSetLogging(solver, 1);

      /* Setup and solve */
      HYPRE_StructSMGSetup(solver, A, b, x);
      HYPRE_StructSMGSolve(solver, A, b, x);

      /* Get some info on the run */
      HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
      HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

      /* Clean up */
      HYPRE_StructSMGDestroy(solver);
   }
   
   if (solver_id == 2)
   {
      HYPRE_StructJacobiCreate(MPI_COMM_WORLD, &solver);
      HYPRE_StructJacobiSetMaxIter(solver, 2000);
      HYPRE_StructJacobiSetTol(solver, 1.0e-12);
      HYPRE_StructJacobiSetZeroGuess(solver);
      HYPRE_StructJacobiSetup(solver, A, b, x);
      HYPRE_StructJacobiSolve(solver, A, b, x);
      HYPRE_StructJacobiGetNumIterations(solver, &num_iterations);
      HYPRE_StructJacobiGetFinalRelativeResidualNorm(solver, &final_res_norm);
      HYPRE_StructJacobiDestroy(solver);
   }


   /* Save the solution for GLVis visualization, see vis/glvis-ex3.sh */
   if (vis)
   {
      FILE *file;
      char filename[255];

      int nvalues = nx*ny;
      double *values = calloc(nvalues, sizeof(double));
      double *exact_solution = calloc(nvalues, sizeof(double));

      /* get the local solution */ 
      HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);
      HYPRE_StructVectorGetBoxValues(analytical, ilower, iupper, exact_solution);
      double diff, diff2,L2;
      double sum = 0;
      for(i=0;i<nvalues;i++){
         diff = values[i] - exact_solution[i];
         diff2 = diff * diff;
         sum += diff2;
         }
      L2 = sqrt(sum/nvalues);
      printf("The L2 Norm is %e\n", L2);
      sprintf(filename, "%s.%06d", "vis/try2varxy.sol", myid);
      if ((file = fopen(filename, "w")) == NULL)
      {
         printf("Error: can't open output file %s\n", filename);
         MPI_Finalize();
         exit(1);
      }

      /* save solution with global unknown numbers */ 
      k = 0;
      for (j = 0; j < ny; j++)
         for (i = 0; i < nx; i++)
            fprintf(file, "%06d %.14e\n", pj*Ny*ny*nx+pi*nx+j*Ny*ny+i, values[k++]);

      fflush(file);
      fclose(file);
      free(values);

     /*  save global finite element mesh */
      if (myid == 0)
	GLVis_PrintGlobalMesh( "vis/try2varxy.mesh", Nx, Ny, nx, ny, hx, hy);
        
   }

   if (myid == 0)
   {
      printf("\n");
      printf("Iterations = %d\n", num_iterations);
      printf("Final Relative Residual Norm = %g\n", final_res_norm);
      printf("\n");
   }

   HYPRE_StructVectorPrint("try2varxy.vect.xf", x, 0); 

   /* Free memory */
   HYPRE_StructGridDestroy(grid);
   HYPRE_StructStencilDestroy(stencil);
   HYPRE_StructMatrixDestroy(A);
   HYPRE_StructVectorDestroy(b);
   HYPRE_StructVectorDestroy(analytical);
   HYPRE_StructVectorDestroy(x);

   /* Finalize MPI */
   MPI_Finalize();

   return (0);
}
