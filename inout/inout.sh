#!/bin/bash
#
gfortran -c -g inout.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling inout.f"
  exit
fi
rm compiler.txt
#
gfortran inout.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading inout.o"
  exit
fi
rm inout.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/inout
#
echo "Executable installed as ~/binf77/$ARCH/inout"
