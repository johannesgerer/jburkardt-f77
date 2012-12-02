#!/bin/bash
#
gfortran -c -g bisection_integer_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bisection_integer_prb.f"
  exit
fi
rm compiler.txt
#
gfortran bisection_integer_prb.o -L$HOME/libf77/$ARCH -lbisection_integer
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bisection_integer_prb.o"
  exit
fi
rm bisection_integer_prb.o
#
mv a.out bisection_integer_prb
./bisection_integer_prb > bisection_integer_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bisection_integer_prb"
  exit
fi
rm bisection_integer_prb
#
echo "Program output written to bisection_integer_prb_output.txt"
