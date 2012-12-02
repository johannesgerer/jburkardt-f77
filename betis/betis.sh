#!/bin/bash
#
gfortran -c -g betis.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling betis.f"
  exit
fi
rm compiler.txt
#
gfortran betis.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading betis.o"
  exit
fi
rm betis.o
#
mv a.out ~/binf77/$ARCH/betis
#
echo "Executable installed as ~/binf77/$ARCH/betis."
