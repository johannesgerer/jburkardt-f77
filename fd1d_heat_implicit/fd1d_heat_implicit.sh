#!/bin/bash
#
F77 -c -g fd1d_heat_implicit.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_heat_implicit.f90"
  exit
fi
rm compiler.txt
#
F77 fd1d_heat_implicit.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd1d_heat_implicit.o"
  exit
fi
rm fd1d_heat_implicit.o
#
mv a.out ~/binf77/$ARCH/fd1d_heat_implicit
#
echo "Executable installed as ~/binf77/$ARCH/fd1d_heat_implicit"
