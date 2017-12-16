#!/bin/sh
source load.sh
rm -r CMakeCache.txt CMakeFiles/
export TRILINOS_DIR=/cm/shared/apps/trilinos/gcc-4.8.1/openmpi-2.0.1/hdf5-1.18.17/2.10.1/
cmake -DCMAKE_CXX_COMPILER=mpic++ \
-DFFTW_LIB=/cm/shared/apps/fftw/openmpi/gcc/64/3.3.3/lib/libfftw3.a \
-DCMAKE_C_COMPILER=mpicc  ./

