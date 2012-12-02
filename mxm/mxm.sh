#!/bin/bash
#
gfortran -c mxm.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxm.f"
  exit
fi
rm compiler.txt
#
gfortran mxm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxm.o"
  exit
fi
rm mxm.o
#
mv a.out ~/binf77/$ARCH/mxm
#
echo "Executable installed as ~/binf77/$ARCH/mxm"
