#!/bin/bash
#
gfortran -c -g clean77.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling clean77.f"
  exit
fi
rm compiler.txt
#
gfortran clean77.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading clean77.o"
  exit
fi
rm clean77.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/clean77
#
echo "Executable installed as ~/binf77/$ARCH/clean77"
