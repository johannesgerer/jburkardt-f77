#!/bin/bash
#
gfortran -c -g complex_numbers.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling complex_numbers.f"
  exit
fi
rm compiler.txt
#
gfortran complex_numbers.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading complex_numbers.o"
  exit
fi
rm complex_numbers.o
#
chmod ugo+x a.out
mv a.out ~/binf77/$ARCH/complex_numbers
#
echo "Executable installed as ~/binf77/$ARCH/complex_numbers"
