#!/bin/bash
#
gfortran -c -g diaphony.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling diaphony.f"
  exit
fi
rm compiler.txt
#
gfortran diaphony.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading diaphony.o"
  exit
fi
#
rm diaphony.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/diaphony
#
echo "Executable installed as ~/binf77/$ARCH/diaphony"
