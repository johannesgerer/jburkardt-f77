#!/bin/bash
#
gfortran -c -g i4lib_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i4lib_prb.f"
  exit
fi
rm compiler.txt
#
gfortran i4lib_prb.o -L$HOME/libf77/$ARCH -li4lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading i4lib_prb.o"
  exit
fi
rm i4lib_prb.o
#
mv a.out i4lib_prb
./i4lib_prb > i4lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running i4lib_prb"
  exit
fi
rm i4lib_prb
#
echo "Test results written to i4lib_prb_output.txt."
