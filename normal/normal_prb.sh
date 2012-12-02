#!/bin/bash
#
gfortran -c -g normal_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling normal_prb.f"
  exit
fi
rm compiler.txt
#
gfortran normal_prb.o -L$HOME/libf77/$ARCH -lnormal
if [ $? -ne 0 ]; then
  echo "Errors linking and loading normal_prb.o"
  exit
fi
rm normal_prb.o
#
mv a.out normal_prb
./normal_prb > normal_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running normal_prb"
  exit
fi
rm normal_prb
#
echo "Test results written to normal_prb_output.txt."
