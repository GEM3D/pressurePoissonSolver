/*
   Example 8

   Interface:    Semi-Structured interface (SStruct)

   Compile with: make ex8

   Sample run:   mpirun -np 2 ex8

   Description:  This is a two processor example which solves a similar
                 problem to the one in Example 2, and Example 6 (The grid
                 boxes are exactly those in the example diagram in the
                 struct interface chapter of the User's Manual.)

                 The difference with the previous examples is that we use
                 three parts, two with a 5-point and one with a 9-point
                 discretization stencil. The solver is PCG with split-SMG
                 preconditioner.
*/

#include <stdio.h>

/* SStruct linear solvers headers */
#include "HYPRE_sstruct_ls.h"

#include "HYPRE_parcsr_ls.h"
#include "HYPRE_parcsr_mv.h"

#include "vis.c"

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

    int object_type;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    /* Dimensions of each box */
    int nx = 10;
    int ny = 10;
    double h = 1.0/nx;
    double h2 = h*h;

    int ndim = 2;
    int nparts = 2;
    int var = 0;
    int nvars = 1;

    int i;

    int ll[2], lr[2], ul[2], ur[2];
    ll[0] = 1;
    ll[1] = 1;
    lr[0] = nx;
    lr[1] = 1;
    ul[0] = 1;
    ul[1] = ny;
    ur[0] = nx;
    ur[1] = ny;

    /* Ghost indices, needed to set neighbor connectivity */
    int ll_g[2], lr_g[2], ul_g[2], ur_g[2];
    ll_g[0] = ll[0] - 1;
    ll_g[1] = ll[1];
    lr_g[0] = lr[0] + 1;
    lr_g[1] = lr[1];
    ul_g[0] = ul[0] - 1;
    ul_g[1] = ul[1];
    ur_g[0] = ur[0] + 1;
    ur_g[1] = ur[1];

    if (num_procs != 1)
    {
        if (myid ==0) printf("Must run with 1 processor!\n");
        MPI_Finalize();

        return(0);
    }

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
            printf("  -vis : save the solution for GLVis visualization\n");
            printf("\n");
        }

        if (print_usage)
        {
            MPI_Finalize();
            return (0);
        }
    }

    /* -------------------------------------------------------------
       1. Set up the 2D grid.

       This gives the index space in each part. We have one variable
       in each part.
       ------------------------------------------------------------- */

    {
        int part;

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

        for (i = 0; i < nparts; i++)
            HYPRE_SStructGridSetVariables(grid, i, nvars, vartypes);

        /* Describe part connectivity */
        {
            int index_map[2] = {0,1};
            int index_dir[2] = {1,1};
            int part;
            int nbor_part;

            part = 0;
            nbor_part = 1;

            /* Cells just outside of the boundary of part 0 in
               its coordinates */

            HYPRE_SStructGridSetNeighborPart(grid, part, lr_g, ur_g,
                                             nbor_part, ll, ul,
                                             index_map, index_dir);

            part = 1;
            nbor_part = 0;
            HYPRE_SStructGridSetNeighborPart(grid, part, ll_g, ul_g,
                                             nbor_part, lr, ur,
                                             index_map, index_dir);
        }

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

       This determines the non-zero structure of the matrix and
       allows non-stencil relationships between the parts
       ------------------------------------------------------------- */

    {
        int part;

        /* Create the graph object */
        HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, &graph);

        /* See MatrixSetObjectType below */
        object_type = HYPRE_SSTRUCT;
        HYPRE_SStructGraphSetObjectType(graph, object_type);

        /* Use the 5-pt stencil on each part  */
        part = 0;
        HYPRE_SStructGraphSetStencil(graph, part, var, stencil_5pt);

        part = 1;
        HYPRE_SStructGraphSetStencil(graph, part, var, stencil_5pt);

        /* Assemble the graph */
        HYPRE_SStructGraphAssemble(graph);
    }

    /* -------------------------------------------------------------
       4. Set up a SStruct Matrix
       ------------------------------------------------------------- */
    {
        int i,j;
        int part;

        /* Create the empty matrix object */
        HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);

        /* Set the object type (by default HYPRE_SSTRUCT). This determines the
           data structure used to store the matrix.  If you want to use unstructured
           solvers, e.g. BoomerAMG, the object type should be HYPRE_PARCSR.
           If the problem is purely structured (with one part), you may want to use
           HYPRE_STRUCT to access the structured solvers.  Since we have two parts
           with different stencils, we set the object type to HYPRE_SSTRUCT. */
        object_type = HYPRE_SSTRUCT;
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

        {
            /* Fix corner stencils */
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



        /* This is a collective call finalizing the matrix assembly.
           The matrix is now ``ready to be used'' */
        HYPRE_SStructMatrixAssemble(A);

        HYPRE_SStructMatrixPrint("matrix.ex2",A,0);
    }


    /* -------------------------------------------------------------
       5. Set up SStruct Vectors for b and x
       ------------------------------------------------------------- */
    {
        int i;
        int part;

        /* Create an empty vector object */
        HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &x);

        /* As with the matrix,  set the object type for the vectors
           to be the sstruct type */
        object_type = HYPRE_SSTRUCT;
        HYPRE_SStructVectorSetObjectType(b, object_type);
        HYPRE_SStructVectorSetObjectType(x, object_type);

        /* Indicate that the vector coefficients are ready to be set */
        HYPRE_SStructVectorInitialize(b);
        HYPRE_SStructVectorInitialize(x);

        /* Set the vector coefficients over the gridpoints in my first box */
        int nvalues = nx*ny;  /* 6 grid points */
        double *values = calloc(nvalues,sizeof(double));

        {
            /* Set b (RHS) */
            for (i = 0; i < nvalues; i++)
                values[i] = 1.0;

            part = 0;
            HYPRE_SStructVectorSetBoxValues(b, part, ll, ur, var, values);

            part = 1;
            HYPRE_SStructVectorSetBoxValues(b, part, ll, ur, var, values);

            /* Set x (initial guess) */
            for (i = 0; i < nvalues; i++)
                values[i] = 0.0;

            part = 0;
            HYPRE_SStructVectorSetBoxValues(x, part, ll,ur, var, values);

            part = 1;
            HYPRE_SStructVectorSetBoxValues(x, part, ll,ur, var, values);
        }
        free(values);

        /* This is a collective call finalizing the vector assembly.
           The vectors are now ``ready to be used'' */
        HYPRE_SStructVectorAssemble(b);
        HYPRE_SStructVectorAssemble(x);
    }


    /* 6. Set up and use a solver (See the Reference Manual for descriptions
       of all of the options.) */
    {
        /* Create an empty PCG Struct solver */
        HYPRE_SStructPCGCreate(MPI_COMM_WORLD, &solver);

        /* Set PCG parameters */
        HYPRE_SStructPCGSetTol(solver, 1.0e-12 );
        HYPRE_SStructPCGSetPrintLevel(solver, 2);
        HYPRE_SStructPCGSetMaxIter(solver, 500);

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
    }

    /* Save the solution for GLVis visualization, see vis/glvis-ex8.sh */
    if (vis)
    {
        /* Scale indices by h */
        double T[8] = {h,0,0,h,h,0,0,h};

        /* Shift second part by (1,0) */
        double O[4] = {0,0,nx*h,0};

        GLVis_PrintSStructGrid(grid, "vis/ex2.mesh", myid, T, O);
        GLVis_PrintSStructVector(x, 0, "vis/ex2.sol", myid);
        GLVis_PrintData("vis/ex2.data", myid, num_procs);
    }

    /* Free memory */
    HYPRE_SStructGridDestroy(grid);
    HYPRE_SStructStencilDestroy(stencil_5pt);
    HYPRE_SStructGraphDestroy(graph);
    HYPRE_SStructMatrixDestroy(A);
    HYPRE_SStructVectorDestroy(b);
    HYPRE_SStructVectorDestroy(x);

    HYPRE_SStructPCGDestroy(solver);
    HYPRE_SStructSplitDestroy(precond);

    /* Finalize MPI */
    MPI_Finalize();

    return (0);
}
