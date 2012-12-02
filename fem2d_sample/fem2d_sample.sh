#!/bin/bash
#
gfortran -c fem2d_sample.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_sample.f"
  exit
fi
rm compiler.txt
#
gfortran fem2d_sample.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_sample.o"
  exit
fi
#
rm fem2d_sample.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/fem2d_sample
#
echo "Program installed as ~/binf77/$ARCH/fem2d_sample"
