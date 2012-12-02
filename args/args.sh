#!/bin/bash
#
gfortran -c -g args.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling args.f"
  exit
fi
rm compiler.txt
#
gfortran args.o
if [ $? -ne 0 ]; then
  echo "Errors loading args.o"
  exit
fi
rm args.o
#
chmod ugo+x args
mv a.out ~/binf77/$ARCH/args
#
echo "Executable installed as ~/binf77/$ARCH/args."
