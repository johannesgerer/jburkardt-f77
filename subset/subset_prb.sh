#!/bin/bash
#
gfortran -c -g subset_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_prb.f"
  exit
fi
rm compiler.txt
#
gfortran subset_prb.o -L$HOME/libf77/$ARCH -lsubset
if [ $? -ne 0 ]; then
  echo "Errors linking and loading subset_prb.o"
  exit
fi
rm subset_prb.o
#
mv a.out subset_prb
./subset_prb > subset_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running subset_prb"
  exit
fi
rm subset_prb
#
echo "Test results written to subset_prb_output.txt."
