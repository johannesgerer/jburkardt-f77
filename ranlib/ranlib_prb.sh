#!/bin/bash
#
gfortran -c -g ranlib_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ranlib_prb.f"
  exit
fi
rm compiler.txt
#
gfortran ranlib_prb.o -L$HOME/libf77/$ARCH -lranlib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ranlib_prb.o"
  exit
fi
rm ranlib_prb.o
#
mv a.out ranlib_prb
./ranlib_prb > ranlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ranlib_prb"
  exit
fi
rm ranlib_prb
#
echo "Test results written to ranlib_prb_output.txt."
