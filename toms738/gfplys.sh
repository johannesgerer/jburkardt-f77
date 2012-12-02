#!/bin/bash
#
gfortran -c -g gfplys.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error while compiling gfplys.f"
  exit
fi
rm compiler.txt
#
gfortran gfplys.o
if [ $? -ne 0 ]; then
  echo "Error while loading gfplys.o"
  exit
fi
rm gfplys.o
#
mv a.out ~/binf77/$ARCH/gfplys
#
echo "Executable installed as ~/binf77/$ARCH/gflys"
