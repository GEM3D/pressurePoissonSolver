## Compiling and Running

This program uses the [Trilinos Library](https://trilinos.org/).

Once you have the Trilinos library install, you just have to set the 
`TRILINOS_DIR` variable in the bash shell to point to where Trilinos
is installed.

```
$ export TRILNOS_DIR = /trilinos/install/dir/
```

After that is set, it should compile with:

```
$ make all
```

This will create an executable named `heat`.
Trilinos uses MPI, so you will have to run it with `mpirun`:

```
$ mpirun -n [num_procs] ./heat [d_x] [d_y] [n_x] [n_y] 

```
Where:

      d_x                               number of domains in the x direction
      d_y                               number of domains in the y direction
      n_x                               number of cells in the x direction, in
                                        each domain
      n_y                               number of cells in the y direction, in
                                        each domain

With the way that I have it set up now, each domain is contained in it's
own process, so `num_procs` will have to be equal to the number of domains.
