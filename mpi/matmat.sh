#!/bin/bash
#
gfortran -c matmat_mpi.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matmat_mpi.f"
  exit
fi
rm compiler.txt
#
gfortran matmat_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matmat_mpi.o"
  exit
fi
rm matmat_mpi.o
#
mv a.out matmat
mpirun -v -np 4 ./matmat > matmat_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running matmat"
  exit
fi
rm matmat
#
echo "The matmat test problem has been executed."
