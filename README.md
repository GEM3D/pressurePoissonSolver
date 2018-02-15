# pressurePoissonSolver

This site organizes our efforts to develop a pressure Poisson
solver for the GEM3d site.   So far, there are two approaches
we are considering :

* Using HYPRE (LLNL) along with a PCG or AMG solver.

* Using a domain decompostion approach.

# Members of the team :

* Grady Wright
* Donna Calhoun
* Scott Aiton (Undergraduate)
* Brenton Peck (Undergraduate)

# Examples to date :

* hypre/two_part/Gem3d/ex3.c  : Two part example with Neumann boundary conditions.
