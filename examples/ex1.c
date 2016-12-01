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

#include "vis.c"

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
    int nx = 4;
    int ny = 4;
    double h = 1;
    double h2 = h*h;

    int ll[2], lr[2], ul[2], ur[2];
    ll[0] = 1;
    ll[1] = 1;
    lr[0] = nx;
    lr[1] = 1;
    ul[0] = 1;
    ul[1] = ny;
    ur[0] = nx;
    ur[1] = ny;

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

    /* Parse command line */
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

    /* 1. Set up the 2D grid.  This gives the index space in each part.
       We have one variable in each part. */
    {
        int ndim = 2;
        int nparts = 2;
        int part;

        /* Create an empty 2D grid object */
        HYPRE_SStructGridCreate(MPI_COMM_WORLD, ndim, nparts, &grid);

        /* Set the extents of the grid - each processor sets its grid
           boxes.  Each part has its own relative index space numbering. */

        /* Processor 0 owns two boxes - one in part 0 and one in part 1. */
        if (myid == 0)
        {
            /* Add parts */
            {
                part = 0;
                HYPRE_SStructGridSetExtents(grid, part, ll, ur);

                part = 1;
                HYPRE_SStructGridSetExtents(grid, part, ll, ur);
            }
        }

        /* Set the variable type and number of variables on each part. */
        {
            int i;
            int nvars = 1;
            HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};

            for (i = 0; i< nparts; i++)
                HYPRE_SStructGridSetVariables(grid, i, nvars, vartypes);
        }

        /* Now we need to set the spatial relation between each of the parts.
           Since we have the same types of variables on both parts, we can
           use HYPRE_GridSetNeighborPart().  Each processor calls this function
           for each part on which it owns boxes that border a different part. */

        if (myid == 0)
        {
            /* The same cells in part 1's coordinates.  Since we use the same
               index space across all parts, the coordinates coincide. */

            int index_map[2] = {0,1};
            int index_dir[2] = {1,1};
            int part;
            int nbor_part;

            /* Relation between part 0 and part 1 on processor 0 */
            {
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

        /* 2. Define the discretization stencils */
        {
            int ndim = 2;
            int var = 0;
            int entry;

            /* the 5-pt stencil in 2D */
            {
                int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
                int stencil_size = 5;

                HYPRE_SStructStencilCreate(ndim, stencil_size, &stencil_5pt);

                for (entry = 0; entry < 5; entry++)
                    HYPRE_SStructStencilSetEntry(stencil_5pt, entry, offsets[entry], var);
            }
        }

        /* 3. Set up the Graph  - this determines the non-zero structure
           of the matrix and allows non-stencil relationships between the parts */
        {
            int var = 0;
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

        /* 4. Set up a SStruct Matrix */
        {
            int i,j;
            int part;
            int var = 0;

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

            /* Each processor must set the stencil values for their boxes on each part.
               In this example, we only set stencil entries and therefore use
               HYPRE_SStructMatrixSetBoxValues.  If we need to set non-stencil entries,
               we have to use HYPRE_SStructMatrixSetValues. */

            if (myid == 0)
            {
                /* Set the matrix coefficients for some set of stencil entries
                   over all the gridpoints in my first box (account for boundary
                   grid points later) */
                {

                    int nentries = 5;
                    int nvalues  = nx*ny*nentries; /* nx*ny grid points */
                    double *values = calloc(nvalues, sizeof(double));

                    int stencil_indices[5];
                    for (j = 0; j < nentries; j++) /* label the stencil indices -
                                                      these correspond to the offsets
                                                      defined above */
                        stencil_indices[j] = j;

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
            }


            /* For each box, set any coefficients that reach ouside of the
               boundary to 0 */
            if (myid == 0)
            {
                int maxnvalues = nx;  /* Assume that nx== ny */
                double *values = calloc(maxnvalues,sizeof(double));
                int stencil_indices[1];

                for (i = 0; i < maxnvalues; i++)
                    values[i] = 0.0;

                /* Bottom edge of each part */
                stencil_indices[0] = 3;
                part = 0;
                HYPRE_SStructMatrixSetBoxValues(A, part, ll, lr,
                                                var, 1,
                                                stencil_indices, values);
                part = 1;
                HYPRE_SStructMatrixSetBoxValues(A, part, ll, lr,
                                                var, 1,
                                                stencil_indices, values);


                /* Left edge of part 0*/
                stencil_indices[0] = 1;
                part = 0;
                HYPRE_SStructMatrixSetBoxValues(A, part, ll, ul,
                                                var, 1,
                                                stencil_indices, values);


                /* Right edge of part 1 */
                stencil_indices[0] = 2;
                part = 1;
                HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                                var, 1,
                                                stencil_indices, values);

                /* Top edge of parts 0 and 1 */
                stencil_indices[0] = 4;

                part = 0;
                HYPRE_SStructMatrixSetBoxValues(A, part, ul, ur,
                                                var, 1,
                                                stencil_indices, values);

                part = 1;
                HYPRE_SStructMatrixSetBoxValues(A, part, ul, ur,
                                                var, 1,
                                                stencil_indices, values);
                free(values);
            }

            /* This is a collective call finalizing the matrix assembly.
               The matrix is now ``ready to be used'' */
            HYPRE_SStructMatrixAssemble(A);
        }

        HYPRE_SStructMatrixPrint("matrix.ex1",A,1);


        /* 5. Set up SStruct Vectors for b and x */
        {
            int i;
            int part;
            int var = 0;

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

            if (myid == 0)
            {
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
            }

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
            HYPRE_SStructPCGSetMaxIter(solver, 200);

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

            GLVis_PrintSStructGrid(grid, "vis/ex1.mesh", myid, T, O);
            GLVis_PrintSStructVector(x, 0, "vis/ex1.sol", myid);
            GLVis_PrintData("vis/ex1.data", myid, num_procs);
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
}
