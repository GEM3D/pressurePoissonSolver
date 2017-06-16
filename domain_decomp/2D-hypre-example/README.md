## Compiling and Running

This program uses the following third party libraries:

+ fftw3
+ mpi
+ Trilinos

### Configuring
First, you will have to have use `cmake` to generate a makefile.
I added a convient script in order to do this:

```
$ ./configure.sh
```

**Note on Trilinos:**

You have to set the `TRILINOS_DIR` environment variable, to point
to the directory where Trilinos was installed.
In my case it was:

```
$ export TRILINOS_DIR = /home/scott/software/trilinos/
```

### Compiling
After the makefile is generated, it should compile with:

```
$ make all
```

### Running
The executable is named `heat`.
It is used in the following way:

    ./heat {OPTIONS} [d_x] [d_y] [n_x] [n_y]
    
    
      OPTIONS:
    
          -h, --help                        Display this help menu
          d_x                               number of domains in the x direction
          d_y                               number of domains in the y direction
          n_x                               number of cells in the x direction, in
                                            each domain
          n_y                               number of cells in the y direction, in
                                            each domain
          -m[matrix filename]               the file to write the matrix to
          "--" can be used to terminate flag options and force all following
          arguments to be treated as positional options
    
