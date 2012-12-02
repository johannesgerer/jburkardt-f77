#!/bin/bash
#
gfortran -c day1_mpi.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling day1_mpi.f"
  exit
fi
rm compiler.txt
#
gfortran day1_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading day1_mpi.o"
  exit
fi
rm day1_mpi.o
#
mv a.out day1
mpirun -v -np 4 ./day1 > day1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running day1"
  exit
fi
rm day1
#
echo "The day1 test problem has been executed."
