#!/bin/bash
#
gfortran -c -g read_align.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling read_align.f"
  exit
fi
rm compiler.txt
#
gfortran read_align.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading read_align.o"
  exit
fi
rm read_align.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/read_align
#
echo "Program installed as ~/binf77/$ARCH/read_align"
