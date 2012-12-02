#!/bin/bash
#
gfortran -c -g latinize_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling latinize_prb.f"
  exit
fi
rm compiler.txt
#
gfortran latinize_prb.o -L$HOME/libf77/$ARCH -llatinize
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latinize_prb.o"
  exit
fi
rm latinize_prb.o
#
mv a.out latinize_prb
./latinize_prb > latinize_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latinize_prb"
  exit
fi
rm latinize_prb
#
echo "Test program output written to latinize_prb_output.txt."
