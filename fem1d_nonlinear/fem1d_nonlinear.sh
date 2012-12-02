#!/bin/bash
#
gfortran -c -g fem1d_nonlinear.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_nonlinear.f"
  exit
fi
rm compiler.txt
#
gfortran fem1d_nonlinear.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_nonlinear.o"
  exit
fi
rm fem1d_nonlinear.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/fem1d_nonlinear
#
echo "Executable installed as ~/binf77/$ARCH/fem1d_nonlinear"
