#!/bin/bash
#
gfortran -c -g nas.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nas.f"
  exit
fi
rm compiler.txt
#
gfortran nas.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nas.o"
  exit
fi
rm nas.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/nas
#
echo "Executable installed as ~/binf77/$ARCH/nas"
