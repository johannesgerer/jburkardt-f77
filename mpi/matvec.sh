#!/bin/bash
#
gfortran -c matvec_mpi.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matvec_mpi.f"
  exit
fi
rm compiler.txt
#
gfortran matvec_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matvec_mpi.o"
  exit
fi
rm matvec_mpi.o
#
mv a.out matvec
mpirun -v -np 4 ./matvec > matvec_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running matvec"
  exit
fi
rm matvec
#
echo "The matvec test problem has been executed."
