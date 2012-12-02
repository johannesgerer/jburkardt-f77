#!/bin/bash
#
gfortran -c -g fem1d.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d.f"
  exit
fi
rm compiler.txt
#
gfortran fem1d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d.o"
  exit
fi
rm fem1d.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/fem1d
#
echo "Executable installed as ~/binf77/$ARCH/fem1d"
