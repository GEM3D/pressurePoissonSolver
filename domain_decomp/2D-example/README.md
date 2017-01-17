## Compiling and Running

This program uses the following third party libraries:

+ fftw3
+ Eigen3
+ boost

### Configuring
First, you will have to have use `cmake` to generate a makefile.
I added a convient script in order to do this:

```
$ ./configure.sh
```

**Note on Eigen:**

If you downloaded eigen, you will have to tell `cmake` where to look
for it. You can do this by setting the `EIGEN3_ROOT` environment variable,
in my case it was:

```
$ export EIGEN3_ROOT = /home/nvgba/Downloads/eigen-eigen-f562a193118d/
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
          -s[solution filename]             the file to write the solution to
          --cg                              use conjugate gradient for solving gamma
                                            values
          --graph                           use a graph when forming the matrix
          "--" can be used to terminate flag options and force all following
          arguments to be treated as positional options
    
