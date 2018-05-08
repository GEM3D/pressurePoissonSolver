#!/bin/sh
rm -r CMakeCache.txt CMakeFiles/
cmake -DCMAKE_CXX_COMPILER=mpic++ \
-DFFTW_LIB=/cm/shared/apps/fftw/openmpi/gcc/64/3.3.3/lib/libfftw3.a \
-DENABLE_MUELU=true -DENABLE_MUELU_CUDA=true \
-DNUM_GPU=2 \
-DCMAKE_C_COMPILER=mpicc  ./

