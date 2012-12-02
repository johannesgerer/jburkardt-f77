#!/bin/bash
#
gfortran -c -g fem2d_poisson_rectangle_linear.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_poisson_rectangle_linear.f"
  exit
fi
rm compiler.txt
#
gfortran fem2d_poisson_rectangle_linear.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_rectangle_linear.o"
  exit
fi
rm fem2d_poisson_rectangle_linear.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/fem2d_poisson_rectangle_linear
#
echo "Program installed as ~/binf77/$ARCH/fem2d_poisson_rectangle_linear"
