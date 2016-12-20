## Compiling and Running

This program uses the [Eigen library](http://eigen.tuxfamily.org/).
Eigen is a header-only library, so you only have to tell the compiler
where the headers are.

Once you download the Eigen library and extract it, you have to update
the `EigenDir` variable in the `EigenDir.mk` file.

In my case, the `EigenDir.mk` file ended up looking like:
```
EigenDir = /home/nvgba/Downloads/eigen-eigen-f562a193118d/
```

After that is completed, it should compile with:

```
$ make all
```

This will create an executable named `heat`.
It is used in the following way:

```
$ ./heat <num_cells>
```

Where `num_cells` is the number of cells in each direction.
