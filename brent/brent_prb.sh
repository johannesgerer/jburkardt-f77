#!/bin/bash
#
gfortran -c -g brent_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brent_prb.f"
  exit
fi
rm compiler.txt
#
gfortran brent_prb.o -L$HOME/libf77/$ARCH -lbrent
if [ $? -ne 0 ]; then
  echo "Errors linking and loading brent_prb.o"
  exit
fi
rm brent_prb.o
#
mv a.out brent_prb
./brent_prb > brent_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running brent_prb"
  exit
fi
rm brent_prb
#
echo "Test results written to brent_prb_output.txt."
