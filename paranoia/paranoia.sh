#!/bin/bash
#
gfortran -c -g paranoia.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling paranoia.f"
  exit
fi
rm compiler.txt
#
gfortran paranoia.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading paranoia.o"
  exit
fi
rm paranoia.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/paranoia
#
echo "Executable installed as ~/binf77/$ARCH/paranoia"
