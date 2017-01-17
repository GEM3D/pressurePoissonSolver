## Compiling and Running

A makefile is provided. To compile type:

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
