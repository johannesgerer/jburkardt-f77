#!/bin/bash
#
gfortran -c -g test_nonlin_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_nonlin_prb.f"
  exit
fi
rm compiler.txt
#
gfortran test_nonlin_prb.o -L$HOME/libf77/$ARCH -ltest_nonlin
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_nonlin_prb.o"
  exit
fi
rm test_nonlin_prb.o
#
mv a.out test_nonlin_prb
./test_nonlin_prb > test_nonlin_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_nonlin_prb"
  exit
fi
rm test_nonlin_prb
#
echo "Test results written to test_nonlin_prb_output.txt."
