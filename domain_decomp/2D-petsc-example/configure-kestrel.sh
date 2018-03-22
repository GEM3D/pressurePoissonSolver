#!/bin/sh
rm -r CMakeCache.txt CMakeFiles/
cmake -DCMAKE_CXX_COMPILER=mpic++ \
-DFFTW_LIB=/cm/shared/apps/fftw/openmpi/gcc/64/3.3.3/lib/libfftw3.a \
-DENABLE_MUELU=false -DENABLE_MUELU_CUDA=false \
-DZoltan_LIBRARIES_EXTRA="-L/cm/shared/apps/parmetis/4.0.3/lib/ -lparmetis -lmetis  -lm "\
-DNUM_GPU=2 \
-DCMAKE_C_COMPILER=mpicc  ./

