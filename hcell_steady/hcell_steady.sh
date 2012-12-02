#!/bin/bash
#
gfortran -c -g hcell_steady.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hcell_steady.f"
  exit
fi
rm compiler.txt
#
gfortran hcell_steady.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hcell_steady.o"
  exit
fi
rm hcell_steady.o
#
mv a.out ~/binf77/$ARCH/hcell_steady
#
echo "Executable installed as ~/binf77/$ARCH/hcell_steady"
