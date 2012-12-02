#!/bin/bash
#
gfortran -c -g filum_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling filum_prb.f"
  exit
fi
rm compiler.txt
#
gfortran filum_prb.o -L$HOME/libf77/$ARCH -lfilum
if [ $? -ne 0 ]; then
  echo "Errors linking and loading filum_prb.o"
  exit
fi
rm filum_prb.o
#
mv a.out filum_prb
./filum_prb > filum_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running filum_prb"
  exit
fi
rm filum_prb
#
echo "Test results written to filum_prb_output.txt."
