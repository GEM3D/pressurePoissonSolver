## Compiling and Running
This program uses the [Eigen library](http://eigen.tuxfamily.org/).
Eigen is a header-only library, so you only have to tell the compiler
where the headers are.

Once you download the Eigen library and extract it, you have to update
the `EIGEN_DIR` variable in the `EigenDir.mk` file.

In my case, the `EigenDir.mk` file ended up looking like:

```
EIGEN_DIR = /home/nvgba/Downloads/eigen-eigen-f562a193118d/
```

After that is completed, it should compile with:

```
$ make all
```

This will create an executable named `heat`.
It is used in the following way:

```
$ ./heat <num_cells> <num_domains>
```

Where `num_cells` is the number of cells on the entire grid and `num_domains` is the number of
domains to use when solving.
