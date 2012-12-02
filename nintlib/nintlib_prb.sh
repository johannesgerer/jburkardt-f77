#!/bin/bash
#
gfortran -c -g nintlib_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nintlib_prb.f"
  exit
fi
rm compiler.txt
#
gfortran nintlib_prb.o -L$HOME/libf77/$ARCH -lnintlib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nintlib_prb.o"
  exit
fi
rm nintlib_prb.o
#
mv a.out nintlib_prb
./nintlib_prb > nintlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nintlib_prb"
  exit
fi
rm nintlib_prb
#
echo "Test results written to nintlib_prb_output.txt."
