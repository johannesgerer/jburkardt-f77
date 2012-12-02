#!/bin/bash
#
gfortran -c -g genin.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error while compiling genin.f"
  exit
fi
rm compiler.txt
#
gfortran genin.o
if [ $? -ne 0 ]; then
  echo "Error while loading genin.o"
  exit
fi
rm genin.o
#
mv a.out ~/binf77/$ARCH/genin
#
echo "Executable installed as ~/binf77/$ARCH/genin"
