#!/bin/bash
#
gfortran -c -g toeplitz_cholesky_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toeplitz_cholesky_prb.f"
  exit
fi
rm compiler.txt
#
gfortran toeplitz_cholesky_prb.o -L$HOME/libf77/$ARCH -ltoeplitz_cholesky
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toeplitz_cholesky_prb.o"
  exit
fi
rm toeplitz_cholesky_prb.o
#
mv a.out toeplitz_cholesky_prb
./toeplitz_cholesky_prb > toeplitz_cholesky_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toeplitz_cholesky_prb"
  exit
fi
rm toeplitz_cholesky_prb
#
echo "Test program output written to toeplitz_cholesky_prb_output.txt."
