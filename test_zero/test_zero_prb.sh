#!/bin/bash
#
gfortran -c -g test_zero_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_zero_prb.f"
  exit
fi
rm compiler.txt
#
gfortran test_zero_prb.o -L$HOME/libf77/$ARCH -ltest_zero
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_zero_prb.o"
  exit
fi
rm test_zero_prb.o
#
mv a.out test_zero_prb
./test_zero_prb > test_zero_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_zero_prb"
  exit
fi
rm test_zero_prb
#
echo "Test results written to test_zero_prb_output.txt."
