#!/bin/bash
#
gfortran -c -g toms431.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling toms431.f"
  exit
fi
rm compiler.txt
#
gfortran toms431.o
if [ $? -ne 0 ]; then
  echo "Errors while loading toms431.o"
  exit
fi
rm toms431.o
#
mv a.out ~/binf77/$ARCH/toms431
#
echo "Executable installed as ~/binf77/$ARCH/toms431"
