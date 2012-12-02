#!/bin/bash
#
gfortran -c quadrature_mpi.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrature_mpi.f"
  exit
fi
rm compiler.txt
#
gfortran quadrature_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrature_mpi.o"
  exit
fi
rm quadrature_mpi.o
#
mv a.out quadrature
mpirun -v -np 4 ./quadrature > quadrature_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quadrature"
  exit
fi
rm quadrature
#
echo "The quadrature test problem has been executed."
