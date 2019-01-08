# pressurePoissonSolver

This is a library for solving Poisson's equation. 

There are two approaches
* Using BoomerAMG as a preconditioner on Schur compliment Matrix.

* Using a geometric multigrid method.

# Members of the team :

* Grady Wright
* Donna Calhoun
* Scott Aiton

# Required Software
* FFTW
* Zoltan
* PETSc
* CMake

# Compiling
Create a seperate source directory and run cmake in the build directory:
```
cmake /path/to/source
```
Then compile with make:
```
make
```

# Examples:

Both of these examples support both the Schur compliment and geometric multigrid methods:
* apps/2d/steady2d  : Two dimensional example for Poisson's equation
* apps/3d/steady    : Three dimensional example for Poisson's equation
