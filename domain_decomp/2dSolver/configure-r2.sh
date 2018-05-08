#!/bin/sh
export FFTW_DIR=/cm/shared/apps/fftw/intel/3.3.6/
export ZOLTAN_DIR=/home/saiton/zoltan/
export PETSC_DIR=/cm/shared/apps/petsc/11sep17/

rm -r CMakeCache.txt CMakeFiles/
cmake -DCMAKE_CXX_COMPILER=mpic++ \
-DCMAKE_C_COMPILER=mpicc  ./

