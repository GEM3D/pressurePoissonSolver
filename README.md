# thunderegg

This is a library for solving Poisson's equation on adaptively refined block-structured Cartesian grids in two and three dimensions in a distributive parallel environment.

The code offers two solvers:
* Domain decomposition using a Schur complement, which is solved using GMRES preconditioned with BoomerAMG.
 
* Geometric multigrid method using a Fast Adaptive Composite (FAC) type algorithm.

The Schur complement method works well for 2D problems, but not 3D.  The Geometric multigrid method works well for both 2D and 3D.

# Members of the team :

* Scott Aiton
* Donna Calhoun
* Grady Wright

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
