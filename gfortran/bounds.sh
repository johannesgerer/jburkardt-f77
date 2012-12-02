#!/bin/bash
#
gfortran -c -fbounds-check bounds.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bounds.f"
  exit
fi
rm compiler.txt
#
gfortran bounds.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bounds.o"
  exit
fi
rm bounds.o
#
mv a.out bounds
./bounds >& bounds_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bounds."
  exit
fi
rm bounds
#
echo "Program output written to bounds_output.txt"
