#!/bin/bash
#
gfortran -c -g minpack_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling minpack_prb.f"
  exit
fi
rm compiler.txt
#
gfortran minpack_prb.o -L$HOME/libf77/$ARCH -lminpack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading minpack_prb.o"
  exit
fi
rm minpack_prb.o
#
mv a.out minpack_prb
./minpack_prb > minpack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running minpack_prb"
  exit
fi
rm minpack_prb
#
echo "Test results written to minpack_prb_output.txt."
