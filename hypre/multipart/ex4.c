/*
   Example 4

   Interface:    Semi-Structured interface with five parts parts

   Grid structure.
       - Each grid is NxN
       - Use SetNeighborParts for sibling grids (P1-P4)
       - Use GraphAddEntries for connectivity between P0 and P1 and P3.


         |---------|-----|-----|
         |         |     |     |
         |         | P3  | P4  |
         |    P0   |-----|-----|
         |         |     |     |
         |         | P1  | P2  |
         |---------|-----|-----|

*/



#include <stdio.h>

/* SStruct linear solvers headers */
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_krylov.h"

#include "HYPRE_parcsr_ls.h"
#include "HYPRE_parcsr_mv.h"

#include "vis.c"
#include "mpi.h"

#define  PI2 M_PI*M_PI

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a > _b ? _a : _b; })


#include <test_solns.h>

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
    HYPRE_SStructVector   error_vec;
    HYPRE_SStructSolver   solver;
    HYPRE_SStructSolver   precond;

    double **error;
    double **solution;
    solution_t qtrue;
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
      1  : Linear (a*x + b*y)
      2  : Quadratic (a*x^2 + b*y^2)
      3  : Trigonometric (cos(a*2*pi*x)*cos(b*2*pi*y)
      .....
      10 : Original problem
    */
    int solution_type = 1;


    /* Lower Left corners for each part : [ax[n],bx[n]] x [ay[n],by[n]] */
    double ax[5] = {0,1,1.5,1,1.5};
    double ay[5] = {0,0,0,0.5,0.5};

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
    int i, j, n, part;

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
            else if (strcmp(argv[arg_index], "-solution_type") == 0)
            {
                arg_index++;
                solution_type = atoi(argv[arg_index++]);
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
            printf("  -solver_id <N>        : Solver_id (N=0,1,2) [1]\n");
            printf("  -nx <N>               : nx (ny will be set to equal nx) [100]\n");
            printf("  -maxiter <N>          : Maximum number of iterations [500]\n");
            printf("  -tol <t>              : Solver tolerance [1e-10]\n");
            printf("  -print_level <N>      : Solver print level to N (N=0-3) [0]\n");
            printf("  -solution_type <N>    : Solution type (0=constant; 1=linear; 2=quad. 3=trig. [1]\n");
            printf("\n");
        }

        if (print_usage)
        {
            MPI_Finalize();
            return (0);
        }
    }

    if (solver_id < 0 || solver_id > 3)
    {
        printf("Solver_id must be in [0,3]\n");
        MPI_Finalize();
        return(0);
    }

    int ny = nx;
    double h = 1.0/nx;   /* Assume that nx == ny */
    double h2 = h*h;

    double hf = h/2.0;   /* Fine grid mesh size */
    double hf2 = hf*hf;
    double hv[5] = {h,hf,hf,hf,hf};  /* Mesh widths */

    int ndim = 2;       /* Dimension of the problem */
    int nparts = 5;     /* Number of parts (fixed)  */
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
    case 3:
        object_type = HYPRE_SSTRUCT;
        break;
    case 1:
    case 2:
        object_type = HYPRE_PARCSR;    /* Doesn't yet work */
    }

    switch(solution_type)
    {
    case 0:
        qtrue = &qtrue_constant;
        break;
    case 1:
        qtrue = &qtrue_plane;
        break;
    case 2:
        qtrue = &qtrue_parabola;
        break;
    case 3:
        qtrue = &qtrue_trig;
        break;
    case 10:
        /* Use original code */
        break;
    }

    /* -------------------------------------------------------------
       1. Set up the 2D grid.
       ------------------------------------------------------------- */

    {
        /* Create an empty 2D grid object */
        HYPRE_SStructGridCreate(MPI_COMM_WORLD, ndim, nparts, &grid);

        /* Set cell centered data */
        HYPRE_SStructVariable vartypes[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};

        /* Five parts (one coarse; four fine) */
        for(n = 0; n < nparts; n++)
        {
            part = n;
            HYPRE_SStructGridSetExtents(grid, part, ll, ur);
        }

        for(n = 0; n < nparts; n++)
        {
            part = n;
            HYPRE_SStructGridSetVariables(grid, part, nvars, vartypes);
        }

#ifndef USE_GRAPH_ENTRIES
        /* Connectivity between four sibling grids (doesn't yet work)

           Vertical connections   : P1-P2, P3-P4
           Horizontal connections : P1-P3, P2-P4
        */
        {
            int index_dir[2] = {1,1};
            int nbor_part;

            {
                int index_map[2] = {0,1};

                /* P1 and P2 */
                part = 1;
                nbor_part = 2;
                HYPRE_SStructGridSetNeighborPart(grid, part, lr_g, ur_g,
                                                 nbor_part, ll, ul,
                                                 index_map, index_dir);
                part = 2;
                nbor_part = 1;
                HYPRE_SStructGridSetNeighborPart(grid, part, ll_g, ul_g,
                                                 nbor_part, lr, ur,
                                                 index_map, index_dir);
            }
            {
                int index_map[2] = {0,1};

                /* P3 and P4 */
                part = 3;
                nbor_part = 4;
                HYPRE_SStructGridSetNeighborPart(grid, part, lr_g, ur_g,
                                                 nbor_part, ll, ul,
                                                 index_map, index_dir);
                part = 4;
                nbor_part = 3;
                HYPRE_SStructGridSetNeighborPart(grid, part, ll_g, ul_g,
                                                 nbor_part, lr, ur,
                                                 index_map, index_dir);
            }

            {
                int index_map[2] = {1,0};

                /* P1 and P3 */
                part = 1;
                nbor_part = 3;
                HYPRE_SStructGridSetNeighborPart(grid, part, ul_g, ur_g,
                                                 nbor_part, ll, lr,
                                                 index_map, index_dir);
                part = 3;
                nbor_part = 1;
                HYPRE_SStructGridSetNeighborPart(grid, part, ll_g, lr_g,
                                                 nbor_part, ul, ur,
                                                 index_map, index_dir);
            }
            {
                int index_map[2] = {1,0};

                /* P2 and P4 */
                part = 2;
                nbor_part = 4;
                HYPRE_SStructGridSetNeighborPart(grid, part, ul_g, ur_g,
                                                 nbor_part, ll, lr,
                                                 index_map, index_dir);
                part = 4;
                nbor_part = 2;
                HYPRE_SStructGridSetNeighborPart(grid, part, ll_g, lr_g,
                                                 nbor_part, ul, ur,
                                                 index_map, index_dir);
            }
        }
#endif

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
        for(n = 0; n < nparts; n++)
        {
            part = n;
            HYPRE_SStructGraphSetStencil(graph, part, var, stencil_5pt);
        }

#ifdef USE_GRAPH_ENTRIES
        {
            /* Connect sibling grids */
            int idx[2], to_idx[2];
            int to_part;

            {
                /* Vertical connections (P1-P2, P3-P4) */
                {
                    int parts[2] = {1,3};
                    int to_parts[2] = {2,4};

                    for(j = lr[1]; j <= ur[1]; j++)
                    {
                        idx[1] = j;
                        to_idx[1] = j;
                        for(n = 0; n < 2; n++)
                        {
                            part = parts[n];
                            idx[0] = lr[0];

                            to_part = to_parts[n];
                            to_idx[0] = ll[0];

                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);

                            part = to_parts[n];
                            idx[0] = ll[0];

                            to_part = parts[n];
                            to_idx[0] = lr[0];
                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);
                        }
                    }
                }

                /* Horizontal connections (P1-P2, P3-P4) */
                {
                    int parts[2] = {1,2};
                    int to_parts[2] = {3,4};

                    for(i = ul[0]; i <= ur[0]; i++)
                    {
                        idx[0] = i;
                        to_idx[0] = i;
                        for(n = 0; n < 2; n++)
                        {
                            part = parts[n];
                            idx[1] = ul[1];

                            to_part = to_parts[n];
                            to_idx[1] = ll[1];

                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);

                            part = to_parts[n];
                            idx[1] = ll[1];

                            to_part = parts[n];
                            to_idx[1] = ul[1];
                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);
                        }
                    }
                }

                /* Vertical connections (P0-P1 and P0-P3) */
                {
                    /* Lower half */
                    for(j = lr[1]; j <= ny/2; j++)
                    {
                        {
                            /* Lower half */
                            part = 0;
                            idx[0] = lr[0];
                            idx[1] = j;

                            to_part = 1;
                            to_idx[0] = ll[0];
                            to_idx[1] = 2*j;  /* Simple connection - not accurate */

                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);

                            part = 1;
                            idx[0] = ll[0];
                            idx[1] = 2*j;  /* Simple connection - not accurate */

                            to_part = 0;
                            to_idx[0] = lr[0];
                            to_idx[1] = j;
                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);

                            part = 1;
                            idx[0] = ll[0];
                            idx[1] = 2*j-1;  /* Simple connection - not accurate */
                            to_idx[0] = lr[0];
                            to_idx[1] = j;
                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);

                        }
                    }

                    for(j = ny/2+1; j <= ur[1]; j++)
                    {
                        {
                            /* upper half */
                            part = 0;
                            idx[0] = lr[0];
                            idx[1] = j;

                            to_part = 3;
                            to_idx[0] = ll[0];
                            to_idx[1] = 2*j - ny;  /* Simple connection - not accurate */

                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);

                            part = 3;
                            idx[0] = ll[0];
                            idx[1] = 2*j - ny;  /* Simple connection - not accurate */

                            to_part = 0;
                            to_idx[0] = lr[0];
                            to_idx[1] = j;
                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);

                            part = 3;
                            idx[0] = ll[0];
                            idx[1] = 2*j - 1 - ny;  /* Simple connection - not accurate */
                            to_idx[0] = lr[0];
                            to_idx[1] = j;
                            HYPRE_SStructGraphAddEntries(graph, part, idx,
                                                         var, to_part, to_idx, var);

                        }
                    }
                }
            }
        }
#endif

        HYPRE_SStructGraphAssemble(graph);
    }


    /* -------------------------------------------------------------
       4. Set up a SStruct Matrix
       ------------------------------------------------------------- */
    {
        /* Create the empty matrix object */
        HYPRE_SStructMatrixCreate(MPI_COMM_WORLD, graph, &A);

        HYPRE_SStructMatrixSetObjectType(A, object_type);

        /* Get ready to set values */
        HYPRE_SStructMatrixInitialize(A);

        /* Set the matrix coefficients */
        {
            int nentries = 5;
            int nvalues  = nx*ny*nentries; /* nx*ny grid points */
            double *values = calloc(nvalues, sizeof(double));
            double *valuesf = calloc(nvalues, sizeof(double));

            int stencil_indices[5];
            for (j = 0; j < nentries; j++)
            {
                stencil_indices[j] = j;
            }

            /* Set stencil on coarse grid */
            for (i = 0; i < nvalues; i += nentries)
            {
                values[i] = -4.0/h2;
                valuesf[i] = -4.0/hf2;
                for (j = 1; j < nentries; j++)
                {
                    values[i+j] = 1.0/h2;
                    valuesf[i+j] = 1.0/hf2;
                }
            }

            part = 0;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ur,
                                            var, nentries,
                                            stencil_indices, values);

            /* Set same stencil on all fine grids */
            for(n = 1; n < nparts; n++)
            {
                part = n;
                HYPRE_SStructMatrixSetBoxValues(A, part, ll, ur,
                                                var, nentries,
                                                stencil_indices, valuesf);

            }

            free(values);
            free(valuesf);
        }




        /* Set boundary stencils to use homogenous Dirichlet onditions.
           (u0 + u1)/2 = ub
               --> u0 = 2*ub - u1;
               --> u0 + u2 - 2*u1 = u2 - 3*u1 + 2*ub
        */
        {
            int nentries = 2;
            int maxnvalues = nentries*MAX(nx,ny);
            double *values = calloc(maxnvalues,sizeof(double));
            double *valuesf = calloc(maxnvalues,sizeof(double));
            int stencil_indices[2];

            for (i = 0; i < maxnvalues; i += nentries)
            {
                values[i] = -3.0/h2;
                values[i+1] = 0;
                valuesf[i] = -3.0/hf2;
                valuesf[i+1] = 0;
            }

            /* Left edge - part 0 */
            stencil_indices[0] = 0;
            stencil_indices[1] = 1;

            part = 0;
            HYPRE_SStructMatrixSetBoxValues(A, part, ll, ul,
                                            var, nentries,
                                            stencil_indices, values);

            /* Right edge (fine grids)  - P2, P4 */
            {
                stencil_indices[0] = 0;
                stencil_indices[1] = 2;
                int parts[2] = {2,4};

                for(n = 0; n < 2; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                }
            }


            /* Bottom edges - P0, P1 P2 */
            {
                stencil_indices[0] = 0;
                stencil_indices[1] = 3;
                int parts[3] = {0,1,2};
                double *v[3] = {values, valuesf, valuesf};

                for(n = 0; n < 3; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, ll, lr,
                                                    var, nentries,
                                                    stencil_indices, v[n]);
                }
            }

            /* Top edge - P0, P3, P4 */
            {
                stencil_indices[0] = 0;
                stencil_indices[1] = 4;
                int parts[3] = {0, 3, 4};
                double *v[3] = {values, valuesf, valuesf};

                for(n = 0; n < 3; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, ul, ur,
                                                    var, nentries,
                                                    stencil_indices, v[n]);
                }
            }

            free(values);
            free(valuesf);
        }

        {
            /* Zero out stencils set along all edges, since we are using graph entries */
            int nentries = 1;
            int maxnvalues = nentries*MAX(nx,ny);
            double *values = calloc(maxnvalues,sizeof(double));
            double *valuesf = calloc(maxnvalues,sizeof(double));
            int stencil_indices[1];

            for (i = 0; i < maxnvalues; i++)
            {
                values[i] = 0;
            }

            /* Right edge of part 0 (replaced by graph entry) */
            stencil_indices[0] = 2;
            part = 0;
            HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                            var, nentries,
                                            stencil_indices, values);

            /* Left edges of sibling grids (P1,P2,P3,P4) */
            {
                stencil_indices[0] = 1;
                int parts[4] = {1,2,3,4};

                for(n = 0; n < 4; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, ll, ul,
                                                    var, nentries,
                                                    stencil_indices, values);
                }
            }

            /* Right edges of P1,P3*/
            {
                stencil_indices[0] = 2;
                int parts[2] = {1,3};

                for(n = 0; n < 2; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                                    var, nentries,
                                                    stencil_indices, values);
                }
            }

            /* Top edges of P1, P2 */
            {
                stencil_indices[0] = 4;
                int parts[2] = {1,2};

                for(n = 0; n < 2; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, ul, ur,
                                                    var, nentries,
                                                    stencil_indices, values);
                }
            }

            /* Bottom edges of P3, P4*/
            {
                stencil_indices[0] = 3;
                int parts[2] = {3,4};

                for(n = 0; n < 2; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, ll, lr,
                                                    var, nentries,
                                                    stencil_indices, values);
                }
            }


#ifdef USE_GRAPH_ENTRIES
            /* Set the values of the stencil weights for graph entries, added above.
               Each graph entry is the 6th entry in the stencil.
            */
            for (j = 0; j < maxnvalues; j++)
            {
                values[j] = 1.0/h2 ;
                valuesf[j] = 1.0/hf2;
            }

            /* Each stencil has one non-stencil entry - added as 6th entry */
            stencil_indices[0] = 5;

            /* Right edge of P0 */
            part = 0;
            HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                            var, nentries,
                                            stencil_indices, values);

            /* Left edges of P1,P2,P3,P4 */
            {
                int parts[4] = {1,2,3,4};
                for(n = 0; n < 4; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, ll, ul,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                }
            }

            /* Right edges of P1,P3 */
            {
                int parts[2] = {1,3};
                for(n = 0; n < 2; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, lr, ur,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                }
            }
            /* Top edges of P1,P2 */
            {
                int parts[2] = {1,2};
                for(n = 0; n < 2; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, ul, ur,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                }
            }

            /* Bottom edges of P3,P4 */
            {
                int parts[2] = {3,4};
                for(n = 0; n < 2; n++)
                {
                    part = parts[n];
                    HYPRE_SStructMatrixSetBoxValues(A, part, ll, lr,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                }
            }

            /* Fix corner entries at common corner P1,P2,P3,P4 */
            {
                int nentries = 2;
                int stencil_indices[2] = {5,6};
                double valuesf[2] = {1.0/hf2, 1.0/hf2};
                {
                    /* P1 (ur corner) */
                    part = 1;
                    HYPRE_SStructMatrixSetBoxValues(A, part, ur, ur,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                    HYPRE_SStructMatrixSetBoxValues(A, part, ul, ul,
                                                    var, nentries,
                                                    stencil_indices, valuesf);

                }
                {
                    /* P2 (ul corner) */
                    part = 2;
                    HYPRE_SStructMatrixSetBoxValues(A, part, ul, ul,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                }
                {
                    /* P3 (lr corner) */
                    part = 3;
                    HYPRE_SStructMatrixSetBoxValues(A, part, lr,lr,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                    HYPRE_SStructMatrixSetBoxValues(A, part, ll,ll,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                }
                {
                    /* P4 (lr corner) */
                    part = 4;
                    HYPRE_SStructMatrixSetBoxValues(A, part, ll,ll,
                                                    var, nentries,
                                                    stencil_indices, valuesf);
                }
            }


#endif
            free(values);
            free(valuesf);
        }



        /* Fix corner stencils */
        {
            int nentries = 3;                 /* Number of stencil weights to set */
            double values[3] = {-2.0/h2,0,0};    /* setting only 1 stencil per call */
            double valuesf[3] = {-2.0/hf2,0,0};    /* setting only 1 stencil per call */

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
            part = 2;
            stencil_indices[1] = 2;
            stencil_indices[2] = 3;
            HYPRE_SStructMatrixSetBoxValues(A, part, lr, lr,
                                            var, nentries,
                                            stencil_indices, valuesf);

            /* Upper left corner */
            part = 0;
            stencil_indices[1] = 1;
            stencil_indices[2] = 4;
            HYPRE_SStructMatrixSetBoxValues(A, part, ul, ul,
                                            var, nentries,
                                            stencil_indices, values);

            /* Upper right corner */
            part = 4;
            stencil_indices[1] = 2;
            stencil_indices[2] = 4;
            HYPRE_SStructMatrixSetBoxValues(A, part, ur, ur,
                                            var, nentries,
                                            stencil_indices, valuesf);
        }

        /* Set one stencil to Dirichlet */
        {
            int nentries = 3;                 /* Number of stencil weights to set */
            double values[3] = {-6/h2,0,0};    /* setting only 1 stencil per call */

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
        HYPRE_SStructMatrixPrint("matrix.ex4",A,0);

    }

    /* -------------------------------------------------------------
       5. Set up SStruct Vectors for b and x
       ------------------------------------------------------------- */
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
        double **values;
        values = calloc(nparts,sizeof(double*));  /* For each part */
        solution = calloc(nparts,sizeof(double*));

        for (n = 0; n < nparts; n++)
        {
            values[n] = calloc(nvalues,sizeof(double));
            solution[n] = calloc(nvalues,sizeof(double));
        }

        {
            if (solution_type < 10)
            {
                double x,y;
                double q,qx,qy,qxx,qyy;

                for(n = 0; n < nparts; n++)
                {
                    int k = 0;
                    for (j = 0; j < ny; j++)
                    {
                        y = ay[n] + (j+0.5)*hv[n];
                        for (i = 0; i < nx; i++)
                        {
                            x = ax[n] + (i+0.5)*hv[n];

                            /* Set above */
                            qtrue(x,y,&q,&qx,&qy,&qxx,&qyy);
                            values[n][k] = qxx + qyy;
                            solution[n][k] = q;
                            if (n == 0)
                            {
                                /* Left edge of P0 */
                                if (i == 0 && j > 0)
                                {
                                    values[n][k] += qx/hv[n];
                                }
                                else if (j == 0 && i > 0)
                                {
                                    values[n][k] += qy/hv[n];
                                }
                                if (j == ny-1)
                                {
                                    values[n][k] -= qy/hv[n];
                                }
                            }
                            else if (n == 1)
                            {
                                if (j == 0)
                                {
                                    values[n][k] += qy/hv[n];
                                }
                            }
                            else if (n == 2)
                            {
                                /* Right edge  of P2 */
                                if (j == 0)
                                {
                                    values[n][k] += qy/hv[n];
                                }
                                if (i == nx-1)
                                {
                                    values[n][k] -= qx/hv[n];
                                }
                            }
                            else if (n == 3)
                            {
                                /* Top edge */
                                if (j == ny - 1)
                                {
                                    values[n][k] -= qy/hv[n];
                                }
                            }
                            else if (n == 4)
                            {
                                if (j == ny-1)
                                {
                                    values[n][k] -= qy/hv[n];
                                }
                                if (i == nx-1)
                                {
                                    values[n][k] -= qx/hv[n];
                                }
                            }
                            k++;
                        }
                    }
                }
                /* Add Dirichlet condition at lower left corner of P0 */
                double q_xface, q_yface;

                /* x-face solution */
                x = ax[0];
                y = ay[0] + hv[0]/2.0;
                qtrue(x,y,&q_xface,&qx,&qy,&qxx,&qyy);

                /* y-face solution */
                x = ax[0] + hv[0]/2.0;
                y = ay[0];
                qtrue(x,y,&q_yface,&qx,&qy,&qxx,&qyy);
                values[0][0] -= (2*q_xface + 2*q_yface)/h2;
            }
            else
            {
                double x,y;
                double C = 2;  /* Coefficients used in exact solution */
                double D = 2;

                for(n = 0; n < nparts; n++)
                {
                    int k = 0;
                    for (j = 0; j < ny; j++)
                    {
                        y = ay[n] + (j+0.5)*hv[n];
                        for (i = 0; i < nx; i++)
                        {
                            x = ax[n] + (i+0.5)*hv[n];
                            solution[n][k] = cos(C*M_PI*x)*cos(D*M_PI*y);
                            values[n][k] = -(C*C + D*D)*M_PI*M_PI*solution[n][k];

                            k++;
                        }
                    }
                }
                /* Add Dirichlet condition at boundary edge */
                double bx = cos(D*M_PI*h/2.0);    /* xface */
                double by = cos(C*M_PI*h/2.0);    /* yface */
                values[0][0] -= (2*bx + 2*by)/h2;
            }

            for(n = 0; n < nparts; n++)
            {
                part = n;
                HYPRE_SStructVectorSetBoxValues(b, part, ll, ur, var, values[n]);
                HYPRE_SStructVectorSetBoxValues(x, part, ll, ur, var, solution[n]);
            }
        }

        {
            /* Set x (initial guess) */
            for(n = 0; n < nparts; n++)
            {
                for (i = 0; i < nvalues; i++)
                {
                    values[n][i] = 0.0;
                }
                HYPRE_SStructVectorSetBoxValues(x, part, ll,ur, var, values[n]);
            }
        }
        for(n = 0; n < nparts; n++)
        {
            free(values[n]);
        }
        free(values);

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
        else if (solver_id == 3)
        {

            HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, &solver);
            HYPRE_BiCGSTABSetMaxIter((HYPRE_Solver) solver, maxiter );
            HYPRE_BiCGSTABSetAbsoluteTol((HYPRE_Solver) solver, tol );
            HYPRE_BiCGSTABSetPrintLevel((HYPRE_Solver) solver, print_level );
            HYPRE_BiCGSTABSetLogging((HYPRE_Solver) solver, 1 );

#if 0
            precond = NULL;
            HYPRE_SStructBiCGSTABSetPrecond(solver,
                                            HYPRE_SStructDiagScale,
                                            HYPRE_SStructDiagScaleSetup,
                                            precond);
#endif

#if 0
            /* Create a split SStruct solver for use as a preconditioner */
            HYPRE_SStructSplitCreate(MPI_COMM_WORLD, &precond);
            HYPRE_SStructSplitSetMaxIter(precond, 1);
            HYPRE_SStructSplitSetTol(precond, 0.0);
            HYPRE_SStructSplitSetZeroGuess(precond);

            /* Set the preconditioner type to split-SMG */
            HYPRE_SStructSplitSetStructSolver(precond, HYPRE_Jacobi);

            HYPRE_SStructBiCGSTABSetPrecond(solver, HYPRE_SStructSplitSolve,
                                       HYPRE_SStructSplitSetup, precond);
#endif
#if 0
            /* Set the AMG preconditioner parameters */
            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetStrongThreshold(precond, .25);
            HYPRE_BoomerAMGSetCoarsenType(precond, 6);
            HYPRE_BoomerAMGSetTol(precond, 0.0);
            HYPRE_BoomerAMGSetMaxIter(precond, 1);
            HYPRE_BoomerAMGSetPrintLevel(precond, 0);

            /* Set the preconditioner */
            HYPRE_SStructBiCGSTABSetPrecond(solver,
                                            HYPRE_BoomerAMGSolve,
                                            HYPRE_BoomerAMGSetup,
                                            precond);

#endif

            HYPRE_SStructBiCGSTABSetup(solver, A, b, x );
            HYPRE_SStructBiCGSTABSolve(solver, A, b, x);

            HYPRE_SStructBiCGSTABGetNumIterations(solver, &num_iterations);
            HYPRE_SStructBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);


            HYPRE_SStructBiCGSTABDestroy(solver);

        }
        double error_norm[3] = {0,0,0};
        {
            /* Compute error */
            int nvalues = nx*ny;
            int n,i;
            double *xvalues;

            error = calloc(nparts,sizeof(double*));

            /* get the local solution */
            if (object_type == HYPRE_SSTRUCT)
            {
                xvalues = calloc(nparts*nvalues,sizeof(double));

                for(n = 0; n < nparts; n++)
                {
                    part = n;
                    HYPRE_SStructVectorGetBoxValues(x, part, ll, ur, var, &xvalues[n*nvalues]);
                }
            }
            else if (object_type == HYPRE_PARCSR)
            {
                xvalues = hypre_VectorData(hypre_ParVectorLocalVector(par_x));
            }

            for(n = 0; n < nparts; n++)
            {
                error[n] = calloc(nparts*nvalues,sizeof(double));
            }

            double hv[5] = {h2, hf2, hf2, hf2, hf2};
            for (i = 0; i < nvalues; i++)
            {
                for (n = 0; n < nparts; n++)
                {
                    double err;
                    err = fabs(solution[n][i] - xvalues[nvalues*n + i]);
                    error_norm[0] += err*hv[n];
                    error_norm[1] += err*err*hv[n];
                    error_norm[2] = MAX(err,error_norm[2]);
                    error[n][i] = err;
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
        /* Each part is transformed using T*[x;y] + O */
        double *O = calloc(2*nparts,sizeof(double));
        double *T = calloc(4*nparts,sizeof(double));
        for(n = 0; n < nparts; n++)
        {
            /* Transformation */
            T[4*n + 0] = hv[n];
            T[4*n + 1] = 0;
            T[4*n + 2] = 0;
            T[4*n + 3] = hv[n];

            /* Lower left corners; shifted by h (vis assumes index space is [0,n-1] */
            O[2*n + 0] = ax[n] - hv[n];
            O[2*n + 1] = ay[n] - hv[n];
        }

        GLVis_PrintSStructGrid(grid, "vis/ex4.mesh", myid, T, O);

        if (object_type == HYPRE_SSTRUCT)
        {
            HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, &error_vec);
            HYPRE_SStructVectorSetObjectType(error_vec, object_type);
            HYPRE_SStructVectorInitialize(error_vec);
            for(n = 0; n < nparts; n++)
            {
                part = n;
                HYPRE_SStructVectorSetBoxValues(error_vec, part, ll, ur, var, solution[n]);
            }
            HYPRE_SStructVectorAssemble(error_vec);

            GLVis_PrintSStructVector(x, 0, "vis/ex4.sol", myid);
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

            sprintf(filename, "%s.%06d", "vis/ex4.sol", myid);
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
        free(T);
        free(O);


        GLVis_PrintData("vis/ex4.data", myid, num_procs);
    }

    /* Free memory */
    HYPRE_SStructGridDestroy(grid);
    HYPRE_SStructStencilDestroy(stencil_5pt);
    HYPRE_SStructGraphDestroy(graph);
    HYPRE_SStructMatrixDestroy(A);
    HYPRE_SStructVectorDestroy(b);
    HYPRE_SStructVectorDestroy(x);

    for(n = 0; n < nparts; n++)
    {
        free(error[n]);
        free(solution[n]);
    }
    free(error);
    free(solution);


    /* Finalize MPI */
    MPI_Finalize();

    return (0);
}
