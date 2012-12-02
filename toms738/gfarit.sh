#!/bin/bash
#
gfortran -c -g gfarit.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error while compiling gfarit.f"
  exit
fi
rm compiler.txt
#
gfortran gfarit.o
if [ $? -ne 0 ]; then
  echo "Error while loading gfarit.o"
  exit
fi
rm gfarit.o
#
mv a.out ~/binf77/$ARCH/gfarit
#
echo "Executable installed as ~/binf77/$ARCH/gfarit"
