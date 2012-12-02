#!/bin/bash
#
gfortran -c -g spread.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spread.f"
  exit
fi
rm compiler.txt
#
gfortran spread.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spread.o"
  exit
fi
rm spread.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/spread
#
echo "Program installed as ~/binf77/$ARCH/spread"
