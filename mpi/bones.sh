#!/bin/bash
#
gfortran -c bones_mpi.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bones_mpi.f"
  exit
fi
rm compiler.txt
#
gfortran bones_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bones_mpi.o"
  exit
fi
rm bones_mpi.o
#
mv a.out bones
mpirun -v -np 4 ./bones > bones_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bones"
  exit
fi
rm bones
#
echo "The bones test problem has been executed."
