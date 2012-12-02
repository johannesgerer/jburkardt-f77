#!/bin/bash
#
gfortran -c -g hcell.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hcell.f"
  exit
fi
rm compiler.txt
#
gfortran hcell.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hcell.o"
  exit
fi
rm hcell.o
#
mv a.out ~/binf77/$ARCH/hcell
#
echo "Executable installed as ~/binf77/$ARCH/hcell."
