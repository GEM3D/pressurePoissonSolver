/*
   Example 3

   Interface:    Semi-Structured interface with two parts

*/

#include <stdio.h>

/* SStruct linear solvers headers */
#include "HYPRE_sstruct_ls.h"

#include "HYPRE_parcsr_ls.h"
#include "HYPRE_parcsr_mv.h"

#include "vis.c"
#include "mpi.h"

#define  PI2 M_PI*M_PI

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a > _b ? _a : _b; })

int main (int argc, char *argv[])
{
    int myid, num_procs;

    int vis = 0;

    HYPRE_SStructGrid     grid;
    HYPRE_SStructGraph    graph;
    HYPRE_SStructStencil  stencil_5pt;
    HYPRE_SStructMatrix   A;
    HYPRE_SStructVector   b;
    HYPRE_SStructVector   x;
    HYPRE_SStructSolver   solver;
    HYPRE_SStructSolver   precond;

    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;


    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    /* ----------------------------------------------------
       Default values for input parameters
       ---------------------------------------------------- */

    int nx = 200;  /* ny set to nx */

    /*
      0 : PCG (no preconditioning)
      1 : PCG (BoomerAMG preconditioning);
      2 : BoomerAMG (no preconditioning)
    */
    int solver_id = 0;

    int maxiter = 500;    /* Maximum number of iterations */
    double tol = 1e-14;   /* Solver tolerance */

    /*
      Print levels for solvers.
      0   : No printing
      >0  : Print iterations (depends on solver)
    */
    int print_level = 0;


    /*
      0  : Constant solution
      1  : Variable solution
    */
    int solution_type = 1;


    /* -------------------------------------------------------------
       For this code we assume that all parts are local to a single
       processor.
       ------------------------------------------------------------- */

    if (num_procs != 1)
    {
        printf("Must run with 1 processor!\n");
        MPI_Finalize();

        return(0);
    }

    /* Miscellaneous variables */
    int i, j, part;

    /* -------------------------------------------------------------
       0. Parse the command line
       ------------------------------------------------------------- */
    {
        int arg_index = 0;
        int print_usage = 0;

        while (arg_index < argc)
        {
            if ( strcmp(argv[arg_index], "-vis") == 0 )
            {
                arg_index++;
                vis = 1;
            }
            else if ( strcmp(argv[arg_index], "-help") == 0 )
            {
                print_usage = 1;
                break;
            }
            else if (strcmp(argv[arg_index], "-solver_id") == 0)
            {
                arg_index++;
                solver_id = atoi(argv[arg_index++]);
            }
            else if (strcmp(argv[arg_index], "-nx") == 0)
            {
                arg_index++;
                nx = atoi(argv[arg_index++]);
            }
            else if (strcmp(argv[arg_index], "-maxiter") == 0)
            {
                arg_index++;
                maxiter = atoi(argv[arg_index++]);
            }
            else if (strcmp(argv[arg_index], "-tol") == 0)
            {
                arg_index++;
                tol = atof(argv[arg_index++]);
            }
            else if (strcmp(argv[arg_index], "-print_level") == 0)
            {
                arg_index++;
                print_level = atoi(argv[arg_index++]);
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
            printf("  -vis                  : Save the solution for GLVis visualization\n");
            printf("  -solver_id N          : Solver_id (N=0,1,2) [1]\n");
            printf("  -nx N                 : nx (ny will be set to equal nx) [100]\n");
            printf("  -maxiter N            : Maximum number of iterations [500]\n");
            printf("  -tol t                : Solver tolerance [1e-10]\n");
            printf("  -print_level N        : Solver print level to N (N=0-3) [0]\n");
            printf("\n");
        }

        if (print_usage)
        {
            MPI_Finalize();
            return (0);
        }
    }

    if (solver_id < 0 || solver_id > 2)
    {
        printf("Solver_id must be 0, 1, or 2\n");
        MPI_Finalize();
        return(0);
    }

    int ny = nx;
    double h = 1.0/nx;   /* Assume that nx == ny */
    double h2 = h*h;

    int ndim = 2;       /* Dimension of the problem */
    int nparts = 2;     /* Number of parts (fixed)  */
    int var = 0;        /* Variable index */
    int nvars = 1;      /* Number of variables */

    /* Index space for each part */
    int ll[2], lr[2], ul[2], ur[2];
    ll[0] = 1;
    ll[1] = 1;
    lr[0] = nx;
    lr[1] = 1;
    ul[0] = 1;
    ul[1] = ny;
    ur[0] = nx;
    ur[1] = ny;

    printf("\n");
#ifdef USE_GRAPH_ENTRIES
    printf("Using GraphAddEntries to define connectivity\n");
#else
    printf("Using SetNeighborPart to define connectivity\n");
#endif
    printf("\n");

#ifndef USE_GRAPH_ENTRIES
    /* Ghost indices, needed for SetNeighborPart */
    int ll_g[2], lr_g[2], ul_g[2], ur_g[2];
    ll_g[0] = ll[0] - 1;
    ll_g[1] = ll[1];
    lr_g[0] = lr[0] + 1;
    lr_g[1] = lr[1];
    ul_g[0] = ul[0] - 1;
    ul_g[1] = ul[1];
    ur_g[0] = ur[0] + 1;
    ur_g[1] = ur[1];
#endif

    int object_type;
    switch(solver_id)
    {
    case 0:
        object_type = HYPRE_SSTRUCT;
        break;
    case 1:
    case 2:
        object_type = HYPRE_PARCSR;    /* Doesn't yet work */
    }



    /* -------------------------------------------------------------
       1. Set up the 2D grid.
       ------------------------------------------------------------- */

    {
        /* Create an empty 2D grid object */
        HYPRE_SStructGridCreate(MPI_COMM_WORLD, ndim, nparts, &grid);

        /* Set the extents of the grid - each processor sets its grid
           boxes.  Each part has its own relative index space numbering. */

        part = 0;
        HYPRE_SStructGridSetExtents(grid, part, ll, ur);

        part = 1;
        HYPRE_SStructGridSetExtents(grid, part, ll, ur);

        /* Set cell centered data */
        HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};

        part = 0;
        HYPRE_SStructGridSetVariables(grid, part, nvars, vartypes);

        part = 1;
        HYPRE_SStructGridSetVariables(grid, part, nvars, vartypes);

#ifndef USE_GRAPH_ENTRIES
        /* Describe part connectivity */
        {
            int index_map[2] = {0,1};
            int index_dir[2] = {1,1};
            int nbor_part;

            /* connectivity described in terms of both part and neighbor
               coordinates */
            part = 0;
            nbor_part = 1;
            HYPRE_SStructGridSetNeighborPart(grid, part, lr_g, ur_g,
                                             nbor_part, ll, ul,
                                             index_map, index_dir);

            part = 1;
            nbor_part = 0;
            HYPRE_SStructGridSetNeighborPart(grid, part, ll_g, ul_g,
                                             nbor_part, lr, ur,
                                             index_map, index_dir);
        }
#endif

        /* Now the grid is ready to use */
        HYPRE_SStructGridAssemble(grid);
    }


    /* -------------------------------------------
       2. Define the discretization stencils
       ------------------------------------------- */
    {
        int entry;

        int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
        int stencil_size = 5;

        HYPRE_SStructStencilCreate(ndim, stencil_size, &stencil_5pt);

        for (entry = 0; entry < stencil_size; entry++)
            HYPRE_SStructStencilSetEntry(stencil_5pt, entry, offsets[entry], var);

    }


    /* -------------------------------------------------------------
       3. Set up the graph.
       ------------------------------------------------------------- */
    {
        /* Create the graph object */
        HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);

        /* See MatrixSetObjectType below */
        HYPRE_SStructGraphSetObjectType(graph, object_type);

        /* Use the 5-pt stencil on each part  */
        part = 0;
        HYPRE_SStructGraphSetStencil(graph, part, var, stencil_5pt);

        part = 1;
        HYPRE_SStructGraphSetStencil(graph, part, var, stencil_5pt);

#ifdef USE_GRAPH_ENTRIES
        /* Add additional graph entries to define connectivity */
        {
            int idx[2], to_idx[2];
            int to_part;

            for(j = lr[1]; j <= ur[1]; j++)
            {

                idx[1] = j;
                to_idx[1] = j;

                /* Connect part 0 to part 1 */
                part = 0;
                idx[0] = lr[0];

                to_part = 1;
                to_idx[0] = ll[0];

                HYPRE_SStructGraphAddEntries(graph, part, idx, var, to_part, to_idx, var);

                /* Connect part 1 to part 0 */
                part = 1;
                idx[0] = ll[0];

                to_part = 0;
                to_idx[0] = lr[0];
                HYPRE_SStructGraphAddEntries(graph, part, idx, var, to_part, to_idx, var);
            }
        }
#endif


        /* Assemble the graph */
        HYPRE_SStructGraphAssemble(graph);
    }

    /* -------------------------------------------------------------
       4. Set up a SStruct Matrix
       ------------------------------------------------------------- */
    {
        /* Create the empty matrix object */
        HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);

        /* Set the object type (by default HYPRE_SSTRUCT). This determines the
           data structure used to store the matrix.  If you want to use unstructured
           solvers, e.g. BoomerAMG, the object type should be HYPRE_PARCSR.
           If the problem is purely structured (with one part), you may want to use
           HYPRE_STRUCT to access the structured solvers.  Since we have two parts
           with different stencils, we set the object type to HYPRE_SSTRUCT. */
        HYPRE_SStructMatrixSetObjectType(A, object_type);

        /* Get ready to set values */
        HYPRE_SStructMatrixInitialize(A);

        /* Set the matrix coefficients */
        {
            int nentries = 5;
            int nvalues  = nx*ny*nentries; /* nx*ny grid points */
            double *values = calloc(nvalues, sizeof(double));

            int stencil_indices[5];
            for (j = 0; j < nentries; j++)
            {
                stencil_indices[j] = j;
            }

            for (i = 0; i < nvalues; i += nentries)
            {
                values[i] = 4.0/h2;
                for (j = 1; j < nentries; j++)
                {
                    values[i+j] = -1.0/h2;
                }
            }

            part = 0;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ur,
                                            var, nentries,
                                            stencil_indices, values);
            part = 1;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ur,
                                            var, nentries,
                                            stencil_indices, values);

            free(values);
        }

        /* Set boundary stencils to use homogenous Neumann conditions. */
        {
            int nentries = 2;
            int maxnvalues = nentries*MAX(nx,ny);
            double *values = calloc(maxnvalues,sizeof(double));
            int stencil_indices[2];

            for (i = 0; i < maxnvalues; i += nentries)
            {
                values[i] = 3.0/h2;
                values[i+1] = 0;
            }

            /* Left edge - part 0 */
            part = 0;
            stencil_indices[0] = 0;
            stencil_indices[1] = 1;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ul,
                                            var, nentries,
                                            stencil_indices, values);

            /* Right edge  - part 1 */
            part = 1;
            stencil_indices[0] = 0;
            stencil_indices[1] = 2;
            HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                            var, nentries,
                                            stencil_indices, values);

            /* Bottom edge - part 0 */
            part = 0;
            stencil_indices[0] = 0;
            stencil_indices[1] = 3;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, lr,
                                            var, nentries,
                                            stencil_indices, values);

            /* Bottom edge - part 1 */
            part = 1;
            stencil_indices[0] = 0;
            stencil_indices[1] = 3;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, lr,
                                            var, nentries,
                                            stencil_indices, values);

            /* Top edge - part 0 */
            part = 0;
            stencil_indices[0] = 0;
            stencil_indices[1] = 4;
            HYPRE_SStructMatrixSetBoxValues(A, part, ul, ur,
                                            var, nentries,
                                            stencil_indices, values);

            /* Top edge - part 1 */
            part = 1;
            stencil_indices[0] = 0;
            stencil_indices[1] = 4;
            HYPRE_SStructMatrixSetBoxValues(A, part, ul, ur,
                                            var, nentries,
                                            stencil_indices, values);

            free(values);
        }

#ifdef USE_GRAPH_ENTRIES
        {
            /* Zero out stencils set along center line */
            int nentries = 1;
            int maxnvalues = nentries*MAX(nx,ny);
            double *values = calloc(maxnvalues,sizeof(double));
            int stencil_indices[1];

            for (i = 0; i < maxnvalues; i += nentries)
            {
                values[i] = 0;
            }

            /* Right edge of part 0 (replaced by graph entry) */
            stencil_indices[0] = 2;
            part = 0;
            HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                            var, nentries,
                                            stencil_indices, values);

            /* Left edge of part 1 (replaced by graph entry) */
            stencil_indices[0] = 1;
            part = 1;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ul,
                                            var, nentries,
                                            stencil_indices, values);

            /* Set the values of the stencil weights for graph entries, added above.
               Each graph entry is the 6th entry in the stencil.
            */
            for (j = 0; j < maxnvalues; j++)
            {
                values[j] = -1.0/h2 ;
            }

            stencil_indices[0] = 5;  /* Graph entry was sixth entry added? */

            part = 0;
            HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                            var, nentries,
                                            stencil_indices, values);

            stencil_indices[0] = 5;  /* Graph entry was sixth entry added? */
            part = 1;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ul,
                                            var, nentries,
                                            stencil_indices, values);
            free(values);
        }


#endif

        /* Fix corner stencils */
        {
            int nentries = 3;                 /* Number of stencil weights to set */
            double values[3] = {2/h2,0,0};    /* setting only 1 stencil per call */

            int stencil_indices[3];
            stencil_indices[0] = 0;

            /* Lower left corner */
            part = 0;
            stencil_indices[1] = 1;
            stencil_indices[2] = 3;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ll,
                                            var, nentries,
                                            stencil_indices, values);

            /* Lower right corner */
            part = 1;
            stencil_indices[1] = 2;
            stencil_indices[2] = 3;
            HYPRE_SStructMatrixSetBoxValues(A, part, lr, lr,
                                            var, nentries,
                                            stencil_indices, values);

            /* Upper left corner */
            part = 0;
            stencil_indices[1] = 1;
            stencil_indices[2] = 4;
            HYPRE_SStructMatrixSetBoxValues(A, part, ul, ul,
                                            var, nentries,
                                            stencil_indices, values);

            /* Upper right corner */
            part = 1;
            stencil_indices[1] = 2;
            stencil_indices[2] = 4;
            HYPRE_SStructMatrixSetBoxValues(A, part, ur, ur,
                                            var, nentries,
                                            stencil_indices, values);
        }

        /* Set one stencil to Dirichlet */
        {
            int nentries = 3;                 /* Number of stencil weights to set */
            double values[3] = {6/h2,0,0};    /* setting only 1 stencil per call */

            int stencil_indices[3];
            stencil_indices[0] = 0;
            stencil_indices[1] = 1;  /* Bottom edge */
            stencil_indices[2] = 3;  /* Bottom edge */

            /* Lower left corner set to Dirichlet */
            part = 0;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ll,
                                            var, nentries,
                                            stencil_indices, values);
        }


        HYPRE_SStructMatrixAssemble(A);
        HYPRE_SStructMatrixPrint("matrix.ex3",A,0);

    }



    /* -------------------------------------------------------------
       5. Set up SStruct Vectors for b and x
       ------------------------------------------------------------- */
    double *solution[2];
    double *error;
    {
        /* Create an empty vector object */
        HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &x);

        /* As with the matrix,  set the object type for the vectors
           to be the sstruct type */
        HYPRE_SStructVectorSetObjectType(b, object_type);
        HYPRE_SStructVectorSetObjectType(x, object_type);

        /* Indicate that the vector coefficients are ready to be set */
        HYPRE_SStructVectorInitialize(b);
        HYPRE_SStructVectorInitialize(x);

        /* Set the vector coefficients over the gridpoints in my first box */
        int nvalues = nx*ny;  /* 6 grid points */
        double *values[2];  /* For each part */
        int n;
        for (n = 0; n < nparts; n++)
        {
            values[n] = calloc(nvalues,sizeof(double));
            solution[n] = calloc(nvalues,sizeof(double));
        }

        {
            if (solution_type == 0)
            {
                /* Set b (RHS) */
                double c = 2;   /* Value of constant solution */
                for (i = 0; i < nvalues; i++)
                {
                    values[0][i] = 0.0;
                    values[1][i] = 0.0;
                    solution[0][i] = c;
                    solution[1][i] = c;
                }
                values[0][0] = c*4/h2;
            }
            else if (solution_type == 1)
            {
                double x,y;
                double C = 2;  /* Coefficients used in exact solution */
                double D = 1;
                int k = 0;
                for (j = 0; j < ny; j++)
                {
                    y = (j+0.5)*h;
                    for (i = 0; i < nx; i++)
                    {
                        /* Part 0 */
                        x = (i+0.5)*h;
                        values[0][k] = (C*C + D*D)*M_PI*M_PI*cos(C*M_PI*x)*cos(D*M_PI*y);
                        solution[0][k] = cos(C*M_PI*x)*cos(D*M_PI*y);

                        /* Part 1 */
                        x += 1.0;
                        values[1][k] = (C*C + D*D)*M_PI*M_PI*cos(C*M_PI*x)*cos(D*M_PI*y);
                        solution[1][k] = cos(C*M_PI*x)*cos(D*M_PI*y);
                        k++;
                    }
                }
                /* Add Dirichlet condition at boundary edge */
                double bx = cos(C*M_PI*h/2.0);    /* Solution at x-face of lower left corner cell */
                double by = cos(D*M_PI*h/2.0);  /* Solution at y-face of lower left corner cell */
                values[0][0] += (2*bx + 2*by)/h2;
            }

            part = 0;
            HYPRE_SStructVectorSetBoxValues(b, part, ll, ur, var, values[0]);

            part = 1;
            HYPRE_SStructVectorSetBoxValues(b, part, ll, ur, var, values[1]);
        }

        {

            /* Set x (initial guess) */
            for (i = 0; i < nvalues; i++)
            {
                values[0][i] = 0.0;
                values[1][i] = 0.0;
            }

            part = 0;
            HYPRE_SStructVectorSetBoxValues(x, part, ll,ur, var, values[0]);

            part = 1;
            HYPRE_SStructVectorSetBoxValues(x, part, ll,ur, var, values[1]);
        }
        free(values[0]);
        free(values[1]);

        HYPRE_SStructVectorAssemble(b);
        HYPRE_SStructVectorAssemble(x);

    }


    /* -------------------------------------------------------------
       6. Set up the solver and solve problem
       ------------------------------------------------------------- */
    {
        double t0, t1;
        double final_res_norm;
        int num_iterations;

        t0 = MPI_Wtime();
        if (solver_id == 0)
        {
            /* This seems to work */

            /* Create an empty PCG Struct solver */
            HYPRE_SStructPCGCreate(MPI_COMM_WORLD, &solver);

            /* Set PCG parameters */
            HYPRE_SStructPCGSetTol(solver, tol );
            HYPRE_SStructPCGSetPrintLevel(solver, print_level);
            HYPRE_SStructPCGSetMaxIter(solver, maxiter);

            /* Create a split SStruct solver for use as a preconditioner */
            HYPRE_SStructSplitCreate(MPI_COMM_WORLD, &precond);
            HYPRE_SStructSplitSetMaxIter(precond, 1);
            HYPRE_SStructSplitSetTol(precond, 0.0);
            HYPRE_SStructSplitSetZeroGuess(precond);

            /* Set the preconditioner type to split-SMG */
            HYPRE_SStructSplitSetStructSolver(precond, HYPRE_SMG);

            /* Set preconditioner and solve */
            HYPRE_SStructPCGSetPrecond(solver, HYPRE_SStructSplitSolve,
                                       HYPRE_SStructSplitSetup, precond);

            HYPRE_SStructPCGSetup(solver, A, b, x);
            HYPRE_SStructPCGSolve(solver, A, b, x);

            HYPRE_SStructPCGGetNumIterations(solver, &num_iterations);
            HYPRE_SStructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

            HYPRE_SStructPCGDestroy(solver);
            HYPRE_SStructSplitDestroy(precond);
        }
        else if (solver_id == 1)
        {
            HYPRE_Solver solver;
            HYPRE_Solver precond;

            /* Set the PCG solvers parameters */
            HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
            HYPRE_ParCSRPCGSetTol(solver,tol);
            HYPRE_ParCSRPCGSetPrintLevel(solver, print_level);
            HYPRE_ParCSRPCGSetMaxIter(solver, maxiter);

#if 0
            /* Not sure what these do */
            HYPRE_ParCSRPCGSetTwoNorm(solver, 1);
            HYPRE_ParCSRPCGSetLogging(solver, 3);
#endif

            /* Set the AMG preconditioner parameters */
            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
            HYPRE_BoomerAMGSetCoarsenType(precond, 6);
            HYPRE_BoomerAMGSetTol(precond, 0.0);
            HYPRE_BoomerAMGSetMaxIter(precond, 1);
            HYPRE_BoomerAMGSetPrintLevel(precond, 0);

            /* Set the preconditioner */
            HYPRE_ParCSRPCGSetPrecond(solver,
                                      HYPRE_BoomerAMGSolve,
                                      HYPRE_BoomerAMGSetup,
                                      precond);

            /* Get matrix and vectors in the sparse matrix format */
            HYPRE_SStructMatrixGetObject(A, (void **) &parcsr_A);
            HYPRE_SStructVectorGetObject(b, (void **) &par_b);
            HYPRE_SStructVectorGetObject(x, (void **) &par_x);

            /* Solve the system */
            HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
            HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

            HYPRE_ParCSRPCGGetNumIterations(solver, &num_iterations);
            HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);


            HYPRE_BoomerAMGDestroy(precond);
            HYPRE_ParCSRPCGDestroy(solver);
        }
        else if (solver_id == 2)
        {
            HYPRE_Solver solver;

            printf("print_level %d\n",print_level);
            printf("maxiter %d\n",maxiter);
            printf("tol %12.4e\n",tol);

            HYPRE_BoomerAMGCreate(&solver);
            HYPRE_BoomerAMGSetTol(solver, tol);
            HYPRE_BoomerAMGSetPrintLevel(solver, print_level);
            HYPRE_BoomerAMGSetMaxIter(solver, maxiter);


#if 0
            /* Stick with default values */
            HYPRE_BoomerAMGSetStrongThreshold(solver, .15);
            HYPRE_BoomerAMGSetNumSweeps(solver, 30);
            HYPRE_BoomerAMGSetMaxLevels(solver, 20);
            HYPRE_BoomerAMGSetCoarsenType(solver, 6);
            HYPRE_BoomerAMGSetRelaxType(solver, 6);
            HYPRE_BoomerAMGSetInterpType(solver, 14);
            HYPRE_BoomerAMGSetAggInterpType(solver, 2);
#endif

            HYPRE_SStructMatrixGetObject(A, (void **) &parcsr_A);
            HYPRE_SStructVectorGetObject(b, (void **) &par_b);
            HYPRE_SStructVectorGetObject(x, (void **) &par_x);


            HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
            HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

            //Get Run Information
            HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
            HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

            HYPRE_BoomerAMGDestroy(solver);
        }

        double error_norm[3] = {0,0,0};
        {
            /* Compute error */
            int nvalues = nx*ny;
            int n,i;
            double *xvalues;

            /* get the local solution */
            if (object_type == HYPRE_SSTRUCT)
            {
                xvalues = calloc(nparts*nvalues,sizeof(double));

                part = 0;
                HYPRE_SStructVectorGetBoxValues(x, part, ll, lr, var, xvalues);
                part = 1;
                HYPRE_SStructVectorGetBoxValues(x, part, ll, lr, var, &xvalues[nvalues]);
            }
            else if (object_type == HYPRE_PARCSR)
            {
                xvalues = hypre_VectorData(hypre_ParVectorLocalVector(par_x));
            }
            error = calloc(nparts*nvalues,sizeof(double));


            for (i = 0; i < nvalues; i++)
            {
                for (n = 0; n < nparts; n++)
                {
                    double err;
                    err = fabs(solution[n][i] - xvalues[nvalues*n + i]);
                    error_norm[0] += err*h2;
                    error_norm[1] += err*err*h2;
                    error_norm[2] = MAX(err,error_norm[2]);
                    error[nvalues*n + i] = err;
                }
            }
            error_norm[0] = error_norm[0]/2.0;     /* Divide by area of domain */
            error_norm[1] = sqrt(error_norm[1]/2.0);
            if (object_type == HYPRE_SSTRUCT)
            {
                free(xvalues);
            }

        }

        t1 = MPI_Wtime();
        printf("%20s %12.4e\n","Elapsed time",t1-t0);
        printf("%20s %12.4e\n","Residual norm",final_res_norm);
        printf("%20s %12d\n","Number of iterations",num_iterations);
        printf("%20s %12.4e %12.4e %12.4e\n","Error",error_norm[0],
               error_norm[1],error_norm[2]);
        if (num_iterations >= maxiter)
        {
            printf(" (Maximum number of iterations exceeded)\n");
        }
        else
        {
            printf("\n");
        }

    }

    /* Save the solution for GLVis visualization, see vis/glvis-ex8.sh */
    if (vis)
    {
        /* Scale indices by h */
        double T[8] = {h,0,0,h,h,0,0,h};

        /* Shift second part by (1,0) */
        double O[4] = {0,0,nx*h,0};

        GLVis_PrintSStructGrid(grid, "vis/ex3.mesh", myid, T, O);
        if (object_type == HYPRE_SSTRUCT)
        {
            GLVis_PrintSStructVector(x, 0, "vis/ex3.sol", myid);
        }
        else if (object_type == HYPRE_PARCSR)
        {
            /* Needed if solution is a ParVector */
            FILE *file;
            char filename[255];

            int nvalues = nparts*nx*ny;
            double *values;

            /* get the local solution */
            values = hypre_VectorData(hypre_ParVectorLocalVector(par_x));
#if 0
            /* To plot error */
            values = error;
#endif

            sprintf(filename, "%s.%06d", "vis/ex3.sol", myid);
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
            for (i = 0; i < nvalues; i++)
                fprintf(file, "%.14e\n", values[i]);

            fflush(file);
            fclose(file);
        }


        GLVis_PrintData("vis/ex3.data", myid, num_procs);
    }

    /* Free memory */
    HYPRE_SStructGridDestroy(grid);
    HYPRE_SStructStencilDestroy(stencil_5pt);
    HYPRE_SStructGraphDestroy(graph);
    HYPRE_SStructMatrixDestroy(A);
    HYPRE_SStructVectorDestroy(b);
    HYPRE_SStructVectorDestroy(x);

    /* Finalize MPI */
    MPI_Finalize();

    return (0);
}
