#!/bin/bash
#
gfortran -c mxv.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxv.f"
  exit
fi
rm compiler.txt
#
gfortran mxv.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxv.o"
  exit
fi
rm mxv.o
#
mv a.out ~/binf77/$ARCH/mxv
#
echo "Executable installed as ~/binf77/$ARCH/mxv"
