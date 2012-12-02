#!/bin/bash
#
gfortran -c -g scpack_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling scpack_prb.f"
  exit
fi
rm compiler.txt
#
gfortran scpack_prb.o -L$HOME/libf77/$ARCH -lscpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading scpack_prb.o"
  exit
fi
rm scpack_prb.o
#
mv a.out scpack_prb
./scpack_prb > scpack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running scpack_prb"
  exit
fi
rm scpack_prb
#
echo "Test results written to scpack_prb_output.txt."
