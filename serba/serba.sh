#!/bin/bash
#
gfortran -c -g serba.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling serba.f"
  exit
fi
rm compiler.txt
#
gfortran serba.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading serba.o"
  exit
fi
rm serba.o
#
mv a.out ~/binf77/$ARCH/serba
#
echo "Executable installed as ~/binf77/$ARCH/serba."
