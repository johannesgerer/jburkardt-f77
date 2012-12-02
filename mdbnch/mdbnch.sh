#!/bin/bash
#
gfortran -c mdbnch.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mdbnch.f"
  exit
fi
rm compiler.txt
#
gfortran mdbnch.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mdbnch.o"
  exit
fi
rm mdbnch.o
#
mv a.out ~/binf77/$ARCH/mdbnch
#
echo "Program installed as ~/binf77/$ARCH/mdbnch"
