#!/bin/bash
#
gfortran -c -g gfortran_intrinsics_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gfortran_intrinsics_prb.f"
  exit
fi
rm compiler.txt
#
gfortran gfortran_intrinsics_prb.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gfortran_intrinsics_prb.o"
  exit
fi
rm gfortran_intrinsics_prb.o
#
mv a.out gfortran_intrinsics_prb
./gfortran_intrinsics_prb > gfortran_intrinsics_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gfortran_intrinsics_prb"
fi
rm gfortran_intrinsics_prb
#
echo "Program output written to gfortran_intrinsics_prb_output.txt"
