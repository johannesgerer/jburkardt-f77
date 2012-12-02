#!/bin/bash
#
gfortran -c -g g77_intrinsics.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling g77_intrinsics.f"
  exit
fi
rm compiler.txt
#
gfortran g77_intrinsics.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading g77_intrinsics.o"
  exit
fi
rm g77_intrinsics.o
#
mv a.out g77_intrinsics
./g77_intrinsics > g77_intrinsics_output.txt
rm g77_intrinsics
#
echo "Program output written to g77_intrinsics_output.txt"
