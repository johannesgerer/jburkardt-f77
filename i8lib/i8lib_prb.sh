#!/bin/bash
#
gfortran -c -g i8lib_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i8lib_prb.f"
  exit
fi
rm compiler.txt
#
gfortran i8lib_prb.o -L$HOME/libf77/$ARCH -li8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading i8lib_prb.o"
  exit
fi
rm i8lib_prb.o
#
mv a.out i8lib_prb
./i8lib_prb > i8lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running i8lib_prb"
  exit
fi
rm i8lib_prb
#
echo "Test results written to i8lib_prb_output.txt."
