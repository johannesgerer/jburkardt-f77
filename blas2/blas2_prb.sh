#!/bin/bash
#
gfortran -c -g blas2_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas2_prb.f"
  exit
fi
rm compiler.txt
#
gfortran blas2_prb.o -L$HOME/libf77/$ARCH -lblas2
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas2_prb.o"
  exit
fi
rm blas2_prb.o
#
mv a.out blas2_prb
./blas2_prb > blas2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas2_prb"
  exit
fi
rm blas2_prb
#
echo "Test results written to blas2_prb_output.txt."
