#!/bin/bash
#
gfortran -c -g blas1_d_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_d_prb.f"
  exit
fi
rm compiler.txt
#
gfortran blas1_d_prb.o -L$HOME/libf77/$ARCH -lblas1_d
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas1_d_prb.o"
  exit
fi
rm blas1_d_prb.o
#
mv a.out blas1_d_prb
./blas1_d_prb > blas1_d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas1_d_prb"
  exit
fi
rm blas1_d_prb
#
echo "Test results written to blas1_d_prb_output.txt."
