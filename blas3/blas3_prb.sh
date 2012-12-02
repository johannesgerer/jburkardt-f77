#!/bin/bash
#
gfortran -c -g blas3_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas3_prb.f"
  exit
fi
rm compiler.txt
#
gfortran blas3_prb.o -L$HOME/libf77/$ARCH -lblas3
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas3_prb.o"
  exit
fi
rm blas3_prb.o
#
mv a.out blas3_prb
./blas3_prb > blas3_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas3_prb"
  exit
fi
rm blas3_prb
#
echo "Test results written to blas3_prb_output.txt."
