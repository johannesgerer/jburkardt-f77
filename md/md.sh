#!/bin/bash
#
gfortran -c -g md.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling md.f"
  exit
fi
rm compiler.txt
#
gfortran md.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md.o"
  exit
fi
rm md.o
#
mv a.out ~/binf77/$ARCH/md
#
echo "Program installed as ~/binf77/$ARCH/md"
