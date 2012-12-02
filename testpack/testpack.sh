#!/bin/bash
#
gfortran -c -g testpack.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling testpack.f"
  exit
fi
rm compiler.txt
#
gfortran testpack.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading testpack.o"
  exit
fi
rm testpack.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/testpack
#
echo "Executable installed as ~/binf77/$ARCH/testpack"
