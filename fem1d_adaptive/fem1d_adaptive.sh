#!/bin/bash
#
gfortran -c -g fem1d_adaptive.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_adaptive.f"
  exit
fi
rm compiler.txt
#
gfortran fem1d_adaptive.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_adaptive.o"
  exit
fi
rm fem1d_adaptive.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/fem1d_adaptive
#
echo "Executable installed as ~/binf77/$ARCH/fem1d_adaptive"
