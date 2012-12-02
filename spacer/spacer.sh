#!/bin/bash
#
gfortran -c -g spacer.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spacer.f"
  exit
fi
rm compiler.txt
#
gfortran spacer.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spacer.o"
  exit
fi
rm spacer.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/spacer
#
echo "Program installed as ~/binf77/$ARCH/spacer"
