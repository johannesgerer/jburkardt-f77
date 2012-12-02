#!/bin/bash
#
gfortran -c -g fem1d_pmethod.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_pmethod.f"
  exit
fi
rm compiler.txt
#
gfortran fem1d_pmethod.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_pmethod.o"
  exit
fi
rm fem1d_pmethod.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/fem1d_pmethod
#
echo "Executable installed as ~/binf77/$ARCH/fem1d_pmethod"
