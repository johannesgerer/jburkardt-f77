#!/bin/bash
#
gfortran -c sum_million.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sum_million.f"
  exit
fi
rm compiler.txt
#
gfortran sum_million.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sum_million.o"
  exit
fi
rm sum_million.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/sum_million
#
echo "Program installed as ~/binf77/$ARCH/sum_million"
