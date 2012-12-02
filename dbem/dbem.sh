#!/bin/bash
#
gfortran -c -g dbem.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dbem.f"
  exit
fi
rm compiler.txt
#
gfortran dbem.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dbem.o"
  exit
fi
rm dbem.o
#
mv a.out ~/binf77/$ARCH/dbem
#
echo "Executable installed as ~/binf77/$ARCH/dbem."
