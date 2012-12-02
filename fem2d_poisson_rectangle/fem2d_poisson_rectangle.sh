#!/bin/bash
#
gfortran -c -g fem2d_poisson_rectangle.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_poisson_rectangle.f"
  exit
fi
rm compiler.txt
#
gfortran fem2d_poisson_rectangle.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_rectangle.o"
  exit
fi
rm fem2d_poisson_rectangle.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/fem2d_poisson_rectangle
#
echo "Executable installed as ~/binf77/$ARCH/fem2d_poisson_rectangle"
