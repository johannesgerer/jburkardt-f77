#!/bin/bash
#
gfortran -c poisson_mpi.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_mpi.f"
  exit
fi
rm compiler.txt
#
gfortran poisson_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading poisson_mpi.o"
  exit
fi
rm poisson_mpi.o
#
mv a.out poisson_band
mpirun -v -np 4 ./poisson_band > poisson_band_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running poisson_band"
  exit
fi
rm poisson_band
#
echo "The poisson_band test problem has been executed."
